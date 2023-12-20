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
                A[i,j,k] = 1.0*rand(rng,1:4,2,2)
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

