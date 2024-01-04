"""
Four players. Each player:

min { xi ∈ ℝ² } : ∑ⱼ∑ₖ₌ⱼ₊₁ xj' Aᵢⱼₖ xk
"""
function setup(::Val{:four_player_matrix_game}; edge_list=[], seed=2, kwargs...)

    # Seed 2 gives nice difference between Nash (edge_list=[]) and parallel
    # bilevel (edge_list=[(1,2), (3,4)]) settings.

    rng = MersenneTwister(seed)
    
    x1 = Symbolics.variables(:x1, 1:2)
    x2 = Symbolics.variables(:x2, 1:2)
    x3 = Symbolics.variables(:x3, 1:2)
    x4 = Symbolics.variables(:x4, 1:2)
    x = Dict(1=>x1, 2=>x2, 3=>x3, 4=>x4)
    
    qp_net = QPNet(x1,x2,x3,x4)
  
    A = Dict() 
    for i = 1:4
        for j = 1:4
            for k = j+1:4
                A[i,j,k] = rand(rng,1:4,2,2)
            end
        end
    end

    for i = 1:4

        cons = [x[i]; sum(x[i])]
        lb = [0.0, 0.0, 1.0]
        ub = [Inf, Inf, 1.0]
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)

        cost = 0
        for j = 1:4
            for k = j+1:4
                cost += x[j]'*A[i,j,k]*x[k]
            end
        end
        dvars = x[i]

        QPN.add_qp!(qp_net, cost, [con_id], dvars)
    end

    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end


function search_for_game(seed_range)
    # 4080 results in 8 
    all_edges = [(1,2),(2,3),(3,4)]
    edge_list_ps = powerset(all_edges) |> collect

    num_unique_equilibria = map(seed_range) do seed
        x_opts = map(edge_list_ps) do edge_list
            qpn = setup(:four_player_matrix_game; edge_list, seed)
            try
                ret = solve(qpn, fill(0.5, 8))
                ret.x_opt
            catch e
                nothing
            end
        end
        if any(isnothing.(x_opts))
            @info "$seed => 0"
            0
        else
            num_unique = length(x_opts)
            for i in 1:length(edge_list_ps)
                qpn = setup(:four_player_matrix_game; edge_list=edge_list_ps[i], seed)
                for j in 1:length(x_opts)
                    j == i && continue
                    try 
                        ret = solve(qpn, x_opts[j])
                        if ret.x_opt ≈ x_opts[i]
                            num_unique -= 1
                            break
                        end
                    catch e
                        continue
                    end
                end
            end
            @info "$seed => $num_unique"
            num_unique
        end
    end
    ii = argmax(num_unique_equilibria)
    @info "The best seed was $(seed_range[ii]) with $(num_unique_equilibria[ii]) unique equilibria"
end
