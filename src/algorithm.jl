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
function solve(qpn::QPNet, x_init, parent_level_request=Set{Linear}(), relaxable_inds=Set{Int}();
        level=1,
        rng=MersenneTwister())

    x_in = x_init
    indent = "       "^level
    req_status = isempty(parent_level_request) ? "is empty" : "is present"
    @info "$indent Solve call for level $level. Parent level request $req_status."
    @info "$indent relaxable inds are $relaxable_inds"

    if isempty(parent_level_request)
        at_level_request = Set{Linear}()
        at_level_inds = Set{Int}()
        while true
            @info "$indent Calling solve_base at level $level. at_level_request is:"
            for req in at_level_request
                @info "$indent      $req"
            end
            @info "$indent x before solve is $x_in"
            (; x_opt, Sol, identified_request, x_alts) = solve_base!(qpn, x_in, at_level_request, at_level_inds; level, rng)
            @info "$indent x after solve is $x_opt"
            if level > 1
                @info "$indent number of solution pieces is $(length(collect(Sol)))"
            end

            if isempty(identified_request)
                @info "$indent No new requests were identified. Returning."
                return (; x_opt, Sol, x_alts)
            else
                if x_in ≈ x_opt
                    @info "$indent Found some new requests. Going to update at_level_request"
                    union!(at_level_request, identified_request)
                else
                    at_level_request = identified_request
                    x_in = x_opt
                end
                at_level_inds = union!(at_level_inds, level_indices(qpn, level))
            end
        end
    else
        # TODO I think this needs to iterate, unioning with any new identified
        # requests (perform after above opts always)
        r_inds = copy(relaxable_inds)
        while true
            @info "$indent Calling FINAL solve_base at level $level. at_level_request is now parent_level_request:"
            for req in parent_level_request
                @info "$indent      $req"
            end
            @info "$indent x before solve is $x_in"
            (; x_opt, Sol, identified_request, x_alts) = solve_base!(qpn, x_in, parent_level_request, r_inds; request_comes_from_parent=true, level, rng)
            @info "$indent x_opt after solve is $x_opt"
            if level > 1
                @info "$indent number of solution pieces is $(length(collect(Sol)))"
            end
            if !(x_opt ≈ x_in)
                @info "$indent when trying to satisfy parent request, different solution found. This shouldn't happen. Returning."
                return (; x_opt, Sol, x_alts)
            elseif isempty(identified_request)
                @info "$indent Remaining requests do not propagate to lower levels. Made best effort to satisfy requests at this level. Returning."
                return (; x_opt, Sol, x_alts)
            else
                @info "$indent Was able to propagate requests to lower level (from $level to $(level+1)). Going to update request and resolve at this level (sending found request to lower level)."
                union!(parent_level_request, identified_request)
                union!(r_inds, level_indices(qpn, level))
            end
        end
    end
end

function solve_base!(qpn::QPNet, x_init, request, relaxable_inds;
        level=1,
        request_comes_from_parent=false,
        rng=MersenneTwister())
    
    if level == num_levels(qpn)
        start = time()
        qep = gather(qpn, level)
        (; x_opt, Sol) = solve_qep(qep, x_init, request, relaxable_inds; 
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
        return (; x_opt, Sol, identified_request=Set{Linear}(), x_alts=Vector{Float64}[])
    else
        x = copy(x_init)
        fair_objective = fair_obj(qpn, level) # TODO should fair_objective still be used for all shared_var modes?
        qep = gather(qpn, level)
        level_constraint_ids = vcat(([id for qp in values(qep.qps) if id ∈ qp.constraint_indices] for id in keys(qep.constraints))...)
        sub_inds = sub_indices(qpn, level)
        param_inds = param_indices(qpn, level)

        for iters in 1:qpn.options.max_iters
            ret_low = solve(qpn, x, request, relaxable_inds; level=level+1, rng)

            x = ret_low.x_opt
            w = x[param_inds]

            Sol_low = ret_low.Sol
            x_alts = ret_low.x_alts

                
            set_guide!(Sol_low, fair_objective)
            start = time()
            local_xs = []
            non_local_xs = Vector{Float64}[]
            local_solutions = [] #Vector{LocalAVISolutions}()
            local_regions = Poly[]
            identified_request = Set{Linear}()
            
            all_same = true
            low_feasible = false
            current_fair_value = fair_objective(x)
            current_infeasible = !all( x ∈ qep.constraints[i].poly for i in level_constraint_ids)
            
            trying_alt = false
            for x_alt in x_alts
                alt_feas = all( x_alt ∈ qep.constraints[i].poly for i in level_constraint_ids)
                alt_not_worse = fair_objective(x_alt) ≤ current_fair_value + qpn.options.tol

                if alt_not_worse && !alt_feas
                    @infiltrate
                end
                if alt_not_worse && alt_feas # Enforcing feasibility for alternate point -- should check this.
                    x = x_alt
                    empty!(request)
                    trying_alt = true
                    break
                end
            end
            if trying_alt
                continue
            end

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
            for (e, S) in enumerate(distinct(Sol_low))
                sub_count += 1
                S_keep = simplify(S)
                low_feasible |= (x ∈ S_keep)
                try
                    res = solve_qep(qep, x, request, relaxable_inds, S_keep, sub_inds;
                                    qpn.var_indices,
                                    level,
                                    subpiece_index=e,
                                    request_comes_from_parent,
                                    qpn.options.debug,
                                    qpn.options.high_dimension,
                                    qpn.options.make_requests,
                                    qpn.options.shared_variable_mode,
                                    rng)
                    set_guide!(res.Sol, z->(z-x)'*(z-x))
                    new_fair_value = fair_objective(res.x_opt) # caution using fair_value
                    better_value_found = new_fair_value < current_fair_value - qpn.options.tol
                    same_value_found = new_fair_value < current_fair_value + qpn.options.tol
                    high_level_shift = !(res.x_opt[param_inds] ≈ w)
            
                    if high_level_shift
                        push!(non_local_xs, res.x_opt)
                        continue
                    end

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

                        @warn "Just emptied request!"
                        empty!(request) # TODO Should this happen here? Seems messy with top logic


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
                catch err
                    err_count += 1
                    continue
                end
            end

            if sub_count > 0 && err_count == sub_count
                error("All subpieces infeasible")
            end

            if !current_infeasible && !low_feasible && false
                @warn "Shouldn't be here..., why was this logic ever here??"
                res = solve_qep(qep, x, request, S_keep, sub_inds; qpn.options.high_dimension, qpn.var_indices)
                diff = norm(x-res.x_opt)
                qpn.options.debug && println("Diff :", diff)
                x .= res.x_opt
                all_same = false
            end

            identified_request = setdiff(identified_request, request)

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
            local S
            try
                S = (qpn.options.gen_solution_map || level > 1) ? combine(local_regions, local_solutions, level_dim; show_progress=true) : nothing
            catch err
                @infiltrate
            end
            # TODO is it needed to specify which subpieces constituted S, and check
            # consistency in up-network solves?
            return (; x_opt=x, Sol=S, identified_request, x_alts=non_local_xs)
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
        error("No solutions to combine... length solutions: 0, length regions: $(length(regions))")
    elseif length(solutions) == 1
        PolyUnion(collect(first(solutions)))
    else
        complements = map(complement, regions)
        it = 0
        combined = map(zip(solutions, complements)) do (s, rc)
            it += 1
            PolyUnion([collect(s); rc.polys])
        end
        #combined = [[collect(s); rc] for (s, rc) in zip(solutions, complements)]
        IntersectionRoot(combined, length.(complements), level_dim; show_progress)
        #PolyUnion(vcat([collect(s) for s in solutions]...))
    end
end
