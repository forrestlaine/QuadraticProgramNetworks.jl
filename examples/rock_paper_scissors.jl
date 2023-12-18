function setup(::Val{:rock_paper_scissors}; kwargs...)
    
    x1 = QPN.variables(:x1, 1:3)
    x2 = QPN.variables(:x2, 1:3)
    ϵ = QPN.variable(:ϵ)

    A = [0.0 1 -1;
         -1 0 1;
         1 -1 0]
    
    qp_net = QPNet(x1,x2,ϵ)
  
    lb = [fill(0.0, 3); 1]
    ub = [fill(Inf, 3); 1]
    cons = [x1 - [ϵ, 0 ,0]; sum(x1)]
    con_id_1 = QPN.add_constraint!(qp_net, cons, lb, ub)
    
    cons = [x2 - [ϵ, 0, 0]; sum(x2)]
    con_id_2 = QPN.add_constraint!(qp_net, cons, lb, ub)
    
    cost_1 = x1'*A*x2
    dvars_1 = x1
    con_ids_1 = [con_id_1]
    qp_id1 = QPN.add_qp!(qp_net, cost_1, con_ids_1, dvars_1...)
    
    cost_2 = x1'*A'*x2
    dvars_2 = x2
    con_ids_2 = [con_id_2]
    qp_id2 = QPN.add_qp!(qp_net, cost_2, con_ids_2, dvars_2...)

    cost_3 = (x1 - [0.5, 0.25, 0.25])'*(x1 - [0.5,0.25,0.25])
    con_id_3 = QPN.add_constraint!(qp_net, [ϵ,], [0.0], [1.0])
    qp_id3 = QPN.add_qp!(qp_net, cost_3, [con_id_3], ϵ)
    
    edge_list = [(qp_id3,qp_id1), (qp_id3,qp_id2)]
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end

