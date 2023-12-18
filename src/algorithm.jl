function solve_base!(qpn::QPNet, x_init, request, relaxable_inds;
        level=1,
        request_comes_from_parent=false,
        rng=MersenneTwister())

    x = copy(x_init)
    for iters in 1:qpn.options.max_iters
        @info "Iteration $iters at level $level"
        
        if level < num_levels(qpn)
            ret_low = solve(qpn, x, request, relaxable_inds; level=level+1, rng)
            @info "Resuming iteration $iters at level $level"
            S = ret_low.Sol
            x = ret_low.x_opt
        else
            S = Dict()
        end

        players_at_level = qpn.network_depth_map[level] |> collect
        processing_tasks = map(players_at_level) do id
            #Threads.@spawn process_qp(qpn, id, x, S)
            process_qp(qpn, id, x, S)
        end
        results = fetch.(processing_tasks)
        equilibrium = true
        subpiece_assignments = Dict{Int, Poly}()
        for (i,id) in enumerate(players_at_level)
            r = results[i]
            if !r.solution
                equilibrium = false
                if level < num_levels(qpn)
                    for (child_id,subpiece_id) in r.subpiece_assignments
                        # Even if another player has already indicated that at least
                        # one subpiece for this particular child (which current
                        # player also parents) results in non-equilibrium,
                        # we choose to overwrite with subpiece causing
                        # discontentment for current player.
                        subpiece_assignments[child_id] = S[child_id][subpiece_id]
                    end
                end
            else
                S[id] = results[i].S
            end
        end
        if !equilibrium
            x = solve_qep(qpn, players_at_level, x, subpiece_assignments)
            continue
        else
            return (; x_opt=x, Sol=S, identified_request=Set{Linear}(), x_alts=Vector{Float64}[])
        end
    end
    error("Can't find solution")
end
