function solve_base!(qpn::QPNet, x_init, request, relaxable_inds;
        level=1,
        proj_vectors=Vector{Vector{Float64}}(),
        request_comes_from_parent=false,
        rng=MersenneTwister())

    x = copy(x_init)
    level == 1 && qpn.options.debug_visualize && qpn.visualization_function(x)
    if level == 1 && isempty(proj_vectors)
        foreach(i->push!(proj_vectors, randn(rng, length(x))), 1:qpn.options.num_projections)
    end
    for iters in 1:qpn.options.max_iters
        proj_vals = [x'v for v in proj_vectors]
        @debug "Iteration $iters at level $level. $proj_vals"
        
        if level < num_levels(qpn)
            ret_low = solve(qpn, x, request, relaxable_inds; level=level+1, rng, proj_vectors)
            if !ret_low.solved
                return (; solved=false, x_fail=x, x_opt=nothing)
            end
            @debug "Resuming iteration $iters at level $level"
            S = ret_low.Sol
            x = ret_low.x_opt
        else
            S = Dict()
        end

        players_at_level = qpn.network_depth_map[level] |> collect |> sort
        players_at_child_level = union((qpn.network_edges[i] for i in players_at_level)...) |> collect |> sort
        processing_tasks = map(players_at_level) do id
            #Threads.@spawn process_qp(qpn, id, x, S)
            process_qp(qpn, id, x, S; exploration_vertices=qpn.options.exploration_vertices)
        end
        results = fetch.(processing_tasks)
        equilibrium = true
        subpiece_assignments = Dict(i=>S[i][1] for i in players_at_child_level)
        subpiece_ids = Dict(i=>1 for i in players_at_child_level)
        for (i,id) in enumerate(players_at_level)
            r = results[i]
            if r.failed
                return (; solved=false, x_fail=x, x_opt=nothing)
            end
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
                        subpiece_ids[child_id] = subpiece_id
                    end
                end
            else
                S[id] = level ∈ qpn.options.levels_to_remove_subsets ? results[i].S |> remove_subsets : results[i].S
                if !isnothing(S[id]) && length(S[id]) == 0
                    @infiltrate
                end
                !isnothing(S[id]) && @debug "Solution graph for node $i has $(length(S[id])) pieces."
            end
        end
        if !equilibrium
            @debug "Equilibium not satisfied at level $level, attempting to find one using the following subpiece assignment:"
            @debug "$(subpiece_ids)"
            try
                xnew = solve_qep(qpn, players_at_level, x, subpiece_assignments)
                if norm(xnew-x) < 1e-3
                    @infiltrate
                end
                x = xnew
                @debug "Equilibrium found, updating solution estimate."
                qpn.options.debug_visualize && qpn.visualization_function(x)
            catch e
                @infiltrate
                @debug "Solving error when computing equilibrium with subpiece ids: $subpiece_ids. Returning x, although this is a known non-equilibrium."
                return (; solved=false, x_fail=x, x_opt=nothing)
            end
            continue
        else
            return (; solved=true, x_opt=x, Sol=S, identified_request=Set{Linear}(), x_alts=Vector{Float64}[])
        end
    end
    error("Can't find solution")
end
