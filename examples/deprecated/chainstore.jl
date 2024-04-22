function setup(::Val{:chainstore}; num_towns=10, α=1.0, exploration_vertices=10)
    A = -[2.0 5.0;
          0.0 5.0]
    B = -[2.0 1.0;
          0.0 1.0]

    p = QPN.variables(:p, 1:2, 1:num_towns)  
    q = QPN.variables(:q, 1:2, 1:num_towns)
    vars = []
    for t in 1:num_towns
        append!(vars, q[:,t])
        append!(vars, p[:,t])
    end
    qp_net = QPNet(vars...)

    pma = p[:,1]
    q_players = Dict()
    p_players = Dict()

    for t in 1:num_towns

        # Ma+Pa move
        pma = α*p[:,t] + (1.0-α)*pma
        cost = pma'*B*q[:,t]
        cons = [q[:,t]; sum(q[:,t])]
        lb = [0.0, 0.0, 1.0]
        ub = [Inf, Inf, 1.0]
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        vars = q[:,t]
        q_players[t] = QPN.add_qp!(qp_net, cost, [con_id,], vars...)


        cost = sum(p[:,s]'*A*q[:,s] for s in t:num_towns)
        cons = [p[:, t]; sum(p[:,t])]
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub) 
        vars = p[:,t]
        p_players[t] = QPN.add_qp!(qp_net, cost, [con_id], vars...)
    end

    # Chain moves
    edge_list = [[(q_players[t], p_players[t]) for t in 1:num_towns]; [(p_players[t], q_players[t+1]) for t in 1:num_towns-1]]
    
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=exploration_vertices, debug_visualize=false)
    qp_net.default_initialization = fill(0.5, 4*num_towns)
    qp_net
end
