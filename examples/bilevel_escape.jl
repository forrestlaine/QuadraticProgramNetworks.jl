function setup(::Val{:bilevel_escape}; kwargs...)
    x = QPN.variables(:x, 1:2) 
    y = QPN.variables(:y, 1:2) 

    qp_net = QPNet(x,y)

    cons1 = [y[1]+y[2], y[1]-y[2]]
    cons2 = [x[1]+x[2], x[1]-x[2]]
    lb = [-1.0, -1]
    ub = -lb

    con_id1 = QPN.add_constraint!(qp_net, cons1, lb, ub)
    con_id2 = QPN.add_constraint!(qp_net, cons2, lb.-1, ub.+1)
           
    cost = 0.5*(y[1]-x[1])^2 + 0.5*(y[2]-x[2])^2
    qp_id2 = QPN.add_qp!(qp_net, cost, [con_id1,], y)

    cost = y[1]-x[1]
    qp_id1 = PN.add_qp!(qp_net, cost, [con_id2,], x)
   
    edge_list = [(qp_id1, qp_id2)]

    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end
