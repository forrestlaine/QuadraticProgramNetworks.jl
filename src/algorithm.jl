"""
Request logic:

1. Temporarily ignore any parent-level request. Optimize to local equilibrium (without making requests).
2. Make requests of lower-levels, if opportunities are identified for further improvement of at-level equilibrium, and re-optimize.
3. Forward parent-level request to lower levels. re-optimize using resulting low-level solution pieces, if different. 
4. Form solution map, adhering to parent-level request when possible. 

Lots of inefficiencies: 
    -bottom level should never ignore requests from parent
    -shouldn't have to ignore always, maybe I can assume that only passed if already at optimum
"""
function solve(qpn::QPNet, x_init, parent_level_request=nothing;
        level=1,
        rng=MersenneTwister())
    at_level_request = Set{Linear}()
    x_opt = x_init
    while true
        (; x_opt, Sol, identified_request) = solve_base(qpn, x_opt, at_level_request; level, rng)
        if isempty(identified_request)
            break
        else
            at_level_request = identified_request
        end
    end
    if isnothing(parent_level_request)
        (; x_opt, Sol)
    else
        (; x_opt, Sol) = solve_base(qpn, x_opt, parent_level_request; request_from_parent=true, level, rng)
    end
end

function solve_base(qpn::QPNet, x_init, request; 
        level=1,
        request_from_parent=false,
        rng=MersenneTwister())

    @infiltrate
    if level == num_levels(qpn)
        start = time()
        qep = gather(qpn, level)
        (; x_opt, Sol) = solve_qep(qep, x_init, request; 
                                   qpn.var_indices,
                                   level,
                                   qpn.options.debug, 
                                   qpn.options.high_dimension, 
                                   qpn.options.shared_variable_mode, 
                                   gen_sol=(num_levels(qpn)==1) ? qpn.options.gen_solution_map : true,
                                   rng)
        fin = time()
        qpn.options.debug && println("Level ", level, " took ", fin-start, " seconds.")
        qpn.options.debug && display_debug(level, 1, x_opt, nothing, nothing)
        return (; x_opt, Sol, identified_request=Set{Linear}())
    else
        x = copy(x_init)
        fair_objective = fair_obj(qpn, level) # TODO should fair_objective still be used for all shared_var modes?
        qep = gather(qpn, level)
        level_constraint_ids = vcat(([id for qp in values(qep.qps) if id ∈ qp.constraint_indices] for id in keys(qep.constraints))...)
        sub_inds = sub_indices(qpn, level)

        for iters in 1:qpn.options.max_iters
            ret_low = solve(qpn, x, request; level=level+1, rng)
            #x_low = ret_low.x_opt
            x = ret_low.x_opt
            Sol_low = ret_low.Sol
            set_guide!(Sol_low, fair_objective)
            start = time()
            local_xs = []
            local_solutions = [] #Vector{LocalAVISolutions}()
            local_regions = Poly[]
            identified_request = Set{Linear}()
            
            all_same = true
            low_feasible = false
            current_fair_value = fair_objective(x)
            current_infeasible = !all( x ∈ qep.constraints[i].poly for i in level_constraint_ids)
            sub_count = 0
            throw_count = 0
            err_count = 0
            if qpn.options.debug #&& level+1 < num_levels(qpn)
                @info "Level $level information:"
                feas = current_infeasible ? "infeasible" : "feasible"
                @info "     Current guess is $(feas)."
                @info "     Current fair value is $(current_fair_value)."
                @info "     About to reason about potentially $(potential_length(Sol_low)) pieces (maybe many more, see lower-level logs)."
            end    
            local S_keep
            @infiltrate level == 1
            for (e, S) in enumerate(distinct(Sol_low))
                sub_count += 1
                S_keep = simplify(S)
                low_feasible |= (x ∈ S_keep)
                try
                    res = solve_qep(qep, x, request, S_keep, sub_inds;
                                    qpn.var_indices,
                                    level,
                                    subpiece_index=e,
                                    qpn.options.debug,
                                    qpn.options.high_dimension,
                                    qpn.options.shared_variable_mode,
                                    rng)
                    @infiltrate
                    set_guide!(res.Sol, z->(z-x)'*(z-x))
                    new_fair_value = fair_objective(res.x_opt) # caution using fair_value
                    better_value_found = new_fair_value < current_fair_value - qpn.options.tol
                    same_value_found = new_fair_value < current_fair_value + qpn.options.tol

                    if current_infeasible || better_value_found
                        diff = norm(x-res.x_opt)
                        if qpn.options.debug
                            if current_infeasible
                                @info "      Previously infeasible, but using piece, found feasible solution. Breaking."
                            else
                                @info "      Better value found ($new_fair_value vs $current_fair_value)! Breaking."
                            end
                        end
                        x .= res.x_opt
                        all_same = false #TODO should queue all non-solutions?
                        break
                    else
                        @debug "     Checking agreement."
                        current_agrees_with_piece = any(S -> x ∈ S, res.Sol)
                        if current_agrees_with_piece
                            @debug "     Agrees."
                        else
                            @debug "     Disagrees."
                        end
                    end
                    if current_agrees_with_piece || same_value_found
                        # assumption here is that if solution has same value (for
                        # fair objective(is this right for games?)) then valid
                        # piece. Warning: These pieces may then be NON-LOCAL.
                        # Needed for some problems (e.g. pessimistic
                        # committment).
                        #
                        # TODO: I think this logic needs to be changed.
                        #
                        push!(local_xs, res.x_opt)
                        push!(local_solutions, res.Sol)
                        push!(local_regions, S_keep)
                        union!(identified_request, res.identified_request)
                    else # poor-valued neighbor
                        throw_count += 1
                        continue
                    end
                catch e
                    err_count += 1
                    @infiltrate
                    continue
                end
            end

            if sub_count > 0 && err_count == sub_count
                error("All subpieces infeasible")
            end

            if !current_infeasible && !low_feasible
                res = solve_qep(qep, x, S_keep, sub_inds; qpn.options.high_dimension, qpn.var_indices)
                diff = norm(x-res.x_opt)
                qpn.options.debug && println("Diff :", diff)
                x .= res.x_opt
                all_same = false
            end
            fin = time()
            qpn.options.debug && println("Level ", level, " took ", fin-start, " seconds.") 
            qpn.options.debug && display_debug(level, iters, x, sub_count, throw_count)

            if qpn.options.high_dimension 
                if (level == 1 && iters < qpn.options.high_dimension_max_iters) || (level > 1 && !all_same)
                    continue
                end
            else
                if !all_same
                    continue
                end
            end

            level_dim = length(param_indices(qpn, level))
            @info "Combining pieces"
            S = (qpn.options.gen_solution_map || level > 1) ? combine(local_regions, local_solutions, level_dim; show_progress=true) : nothing
            # TODO is it needed to specify which subpieces constituted S, and check
            # consistency in up-network solves?
            return (; x_opt=x, Sol=S, identified_request)
        end
        error("Can't find solution.")
    end
end

"""
Conustructs the solution set
S := ⋃ₚ ⋂ᵢ Zᵢᵖ

Zᵢᵖ ∈ { Rᵢ', Sᵢ }
where Rᵢ' is the set complement of Rᵢ.
"""
function combine(regions, solutions, level_dim; show_progress=true)
    if length(solutions) == 0
        @error "No solutions to combine... length solutions: 0, length regions: $(length(regions))"
    elseif length(solutions) == 1
        @info "Length solutions is 1, collecting."
        PolyUnion(collect(first(solutions)))
    else
        @info "Forming complements"
        complements = map(complement, regions)
        @info "Forming combinations"
        it = 0
        combined = map(zip(solutions, complements)) do (s, rc)
            it += 1
            @info it
            PolyUnion([collect(s); rc.polys])
        end
        #combined = [[collect(s); rc] for (s, rc) in zip(solutions, complements)]
        @info "Forming intersection"
        IntersectionRoot(combined, length.(complements), level_dim; show_progress)
        #PolyUnion(vcat([collect(s) for s in solutions]...))
    end
end
