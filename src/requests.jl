function solve(qpn::QPNet)
    solve(qpn, qpn.default_initialization)
end


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
        proj_vectors=Vector{Vector{Float64}}(),
        rng=MersenneTwister(1))
    return solve_base!(qpn, x_init, parent_level_request, relaxable_inds; level, rng, proj_vectors)
    x_in = x_init
    indent = "       "^level
    req_status = isempty(parent_level_request) ? "is empty" : "is present"
    qpn.options.debug && @info "$indent Solve call for level $level. Parent level request $req_status."
    qpn.options.debug && @info "$indent relaxable inds are $relaxable_inds"

    if isempty(parent_level_request)
        at_level_request = Set{Linear}()
        at_level_inds = Set{Int}()
        while true
            qpn.options.debug && @info "$indent Calling solve_base at level $level. at_level_request is:"
            for req in at_level_request
                qpn.options.debug && @info "$indent      $req"
            end
            qpn.options.debug && @info "$indent x before solve is $x_in"
            ret = solve_base!(qpn, x_in, at_level_request, at_level_inds; level, rng)
            (; x_opt, Sol, identified_request, x_alts) = solve_base!(qpn, x_in, at_level_request, at_level_inds; level, rng)
            qpn.options.debug && @info "$indent x after solve is $x_opt"
            if level > 1
                qpn.options.debug && @info "$indent number of solution pieces is $(length(collect(Sol)))"
            end

            if isempty(identified_request)
                qpn.options.debug && @info "$indent No new requests were identified. Returning."
                return (; x_opt, Sol, x_alts)
            else
                if x_in ≈ x_opt
                    qpn.options.debug && @info "$indent Found some new requests. Going to update at_level_request"
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
            qpn.options.debug && @info "$indent Calling FINAL solve_base at level $level. at_level_request is now parent_level_request:"
            for req in parent_level_request
                qpn.options.debug && @info "$indent      $req"
            end
            qpn.options.debug && @info "$indent x before solve is $x_in"
            (; x_opt, Sol, identified_request, x_alts) = solve_base!(qpn, x_in, parent_level_request, r_inds; request_comes_from_parent=true, level, rng)
            qpn.options.debug && @info "$indent x_opt after solve is $x_opt"
            if level > 1
                qpn.options.debug && @info "$indent number of solution pieces is $(length(collect(Sol)))"
            end
            if !(x_opt ≈ x_in)
                qpn.options.debug && @info "$indent when trying to satisfy parent request, different solution found. This shouldn't happen. Returning."
                return (; x_opt, Sol, x_alts)
            elseif isempty(identified_request)
                qpn.options.debug && @info "$indent Remaining requests do not propagate to lower levels. Made best effort to satisfy requests at this level. Returning."
                return (; x_opt, Sol, x_alts)
            else
                qpn.options.debug && @info "$indent Was able to propagate requests to lower level (from $level to $(level+1)). Going to update request and resolve at this level (sending found request to lower level)."
                union!(parent_level_request, identified_request)
                union!(r_inds, level_indices(qpn, level))
            end
        end
    end
end
