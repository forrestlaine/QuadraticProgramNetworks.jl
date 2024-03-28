    #        w = x[param_inds]

    #        Sol_low = ret_low.Sol
    #        x_alts = ret_low.x_alts

    #        start = time()
    #        local_xs = []
    #        non_local_xs = Vector{Float64}[]
    #        local_solutions = [] #Vector{LocalAVISolutions}()
    #        local_regions = Poly[]
    #        identified_request = Set{Linear}()
    #        
    #        all_same = true
    #        low_feasible = false
    #        current_fair_value = fair_objective(x)
    #        current_infeasible = !all( x ∈ qep.constraints[i].poly for i in level_constraint_ids)
    #        
    #        trying_alt = false
    #        for x_alt in x_alts
    #            alt_feas = all( x_alt ∈ qep.constraints[i].poly for i in level_constraint_ids)
    #            alt_not_worse = fair_objective(x_alt) ≤ current_fair_value + qpn.options.tol

    #            if alt_not_worse && !alt_feas
    #                @infiltrate
    #            end
    #            if alt_not_worse && alt_feas # Enforcing feasibility for alternate point -- should check this.
    #                x = x_alt
    #                empty!(request)
    #                trying_alt = true
    #                break
    #            end
    #        end
    #        if trying_alt
    #            continue
    #        end

    #        if qpn.options.try_hull
    #            PU = PolyUnion(collect(distinct(Sol_low)))
    #            H = convex_hull(PU) |> simplify
    #            low_feasible = x ∈ PU
    #            @infiltrate
    #            try 
    #                res = solve_qep(qpn, level, x, request, relaxable_inds, H, sub_inds;
    #                                qpn.var_indices,
    #                                subpiece_index=0,
    #                                request_comes_from_parent,
    #                                qpn.options.debug,
    #                                qpn.options.high_dimension,
    #                                qpn.options.make_requests,
    #                                qpn.options.shared_variable_mode,
    #                                rng)

    #                @infiltrate
    #                new_fair_value = fair_objective(res.x_opt) # caution using fair_value
    #                better_value_found = new_fair_value < current_fair_value - qpn.options.tol
    #                same_value_found = new_fair_value < current_fair_value + qpn.options.tol
    #                high_level_shift = !(res.x_opt[param_inds] ≈ w)
    #                if !(high_level_shift || current_infeasible || better_value_found)
    #                    @infiltrate
    #                    current_agrees_with_piece = any(S -> x ∈ S, res.Sol)
    #                    if current_agrees_with_piece || same_value_found
    #                        valid = all(S ⊆ PU for S in res.Sol)
    #                        @infiltrate
    #                    end
    #                end
    #            catch e
    #                @infiltrate
    #            end
    #        end

    #        sub_count = 0
    #        throw_count = 0
    #        err_count = 0
    #        if qpn.options.debug #&& level+1 < num_levels(qpn)
    #            @info "Level $level information:"
    #            feas = current_infeasible ? "infeasible" : "feasible"
    #            @info "     Current guess is $(feas)."
    #            @info "     Current fair value is $(current_fair_value)."
    #            @info "     About to reason about potentially $(potential_length(Sol_low)) pieces (maybe many more, see lower-level logs)."
    #        end    
    #        local S_keep
    #        for (e, S) in enumerate(distinct(Sol_low))
    #            sub_count += 1
    #            S_keep = simplify(S)
    #            low_feasible |= (x ∈ S_keep)
    #            try
    #                res = solve_qep(qpn, level, x, request, relaxable_inds, S_keep, sub_inds;
    #                                qpn.var_indices,
    #                                subpiece_index=e,
    #                                request_comes_from_parent,
    #                                qpn.options.debug,
    #                                qpn.options.high_dimension,
    #                                qpn.options.make_requests,
    #                                qpn.options.shared_variable_mode,
    #                                rng)
    #                new_fair_value = fair_objective(res.x_opt) # caution using fair_value
    #                better_value_found = new_fair_value < current_fair_value - qpn.options.tol
    #                same_value_found = new_fair_value < current_fair_value + qpn.options.tol
    #                high_level_shift = !(res.x_opt[param_inds] ≈ w)
    #        
    #                if high_level_shift
    #                    push!(non_local_xs, res.x_opt)
    #                    continue
    #                end

    #                if current_infeasible || better_value_found
    #                    diff = norm(x-res.x_opt)
    #                    if qpn.options.debug
    #                        if current_infeasible
    #                            @info "      Previously infeasible, but using piece, found feasible solution. Breaking."
    #                        else
    #                            @info "      Better value found ($new_fair_value vs $current_fair_value)! Breaking."
    #                        end
    #                    end
    #                    x .= res.x_opt

    #                    #@warn "Just emptied request!"
    #                    empty!(request) # TODO Should this happen here? Seems messy with top logic


    #                    all_same = false #TODO should queue all non-solutions?
    #                    break
    #                else
    #                    @debug "     Checking agreement."
    #                    current_agrees_with_piece = any(S -> x ∈ S, res.Sol)
    #                    if current_agrees_with_piece
    #                        @debug "     Agrees."
    #                    else
    #                        @debug "     Disagrees."
    #                    end
    #                end
    #                if current_agrees_with_piece || same_value_found
    #                    # assumption here is that if solution has same value (for
    #                    # fair objective(is this right for games?)) then valid
    #                    # piece. Warning: These pieces may then be NON-LOCAL.
    #                    # Needed for some problems (e.g. pessimistic
    #                    # committment).
    #                    #
    #                    # TODO: I think this logic needs to be changed.
    #                    #
    #                    push!(local_xs, res.x_opt)
    #                    push!(local_solutions, res.Sol)
    #                    push!(local_regions, S_keep)
    #                    union!(identified_request, res.identified_request)
    #                else # poor-valued neighbor
    #                    throw_count += 1
    #                    continue
    #                end
    #            catch err
    #                err_count += 1
    #                continue
    #            end
    #        end

    #        if sub_count > 0 && err_count == sub_count
    #            error("All subpieces infeasible")
    #        end

    #        if !current_infeasible && !low_feasible && false
    #            @warn "Shouldn't be here..., why was this logic ever here??"
    #            res = solve_qep(qep, x, request, S_keep, sub_inds; qpn.options.high_dimension, qpn.var_indices)
    #            diff = norm(x-res.x_opt)
    #            qpn.options.debug && println("Diff :", diff)
    #            x .= res.x_opt
    #            all_same = false
    #        end

    #        identified_request = setdiff(identified_request, request)

    #        fin = time()
    #        qpn.options.debug && println("Level ", level, " took ", fin-start, " seconds.") 
    #        qpn.options.debug && display_debug(level, iters, x, sub_count, throw_count)

    #        if qpn.options.high_dimension 
    #            if (level == 1 && iters < qpn.options.high_dimension_max_iters) || (level > 1 && !all_same)
    #                continue
    #            end
    #        else
    #            if !all_same
    #                continue
    #            end
    #        end
    #        @infiltrate 
    #        level_dim = length(param_indices(qpn, level))
    #        local S
    #        try
    #            S = (qpn.options.gen_solution_map || level > 1) ? combine(local_regions, local_solutions, level_dim; show_progress=true) : nothing
    #        catch err
    #            @infiltrate
    #        end
    #        # TODO is it needed to specify which subpieces constituted S, and check
    #        # consistency in up-network solves?
    #        return (; x_opt=x, Sol=S, identified_request, x_alts=non_local_xs)
    #    end
    #    error("Can't find solution.")
    #end
#end
