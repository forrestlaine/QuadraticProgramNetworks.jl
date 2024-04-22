function setup(::Val{:trilevel_escape}; kwargs...)
    x = QPN.variable(:x) 
    y = QPN.variable(:y) 
    z = QPN.variables(:z, 1:2)

    qp_net = QPNet(x,y,z)

    cons_x = [x,]
    lb_x = [-2.0,]
    ub_x = [2.0,]
    con_id_x = QPN.add_constraint!(qp_net, cons_x, lb_x, ub_x)

    cons_y = [y-z[1], y-z[2]]
    lb_y = [0.0, 0]
    ub_y = [Inf, Inf]
    con_id_y = QPN.add_constraint!(qp_net, cons_y, lb_y, ub_y)

    cons_z1 = [z[1],]
    lb_z1 = [-1.0,]
    ub_z1 = [1.0,]
    con_id_z1 = QPN.add_constraint!(qp_net, cons_z1, lb_z1, ub_z1)

    cons_z2 = [z[2],]
    lb_z2 = [-0.5,]
    ub_z2 = [0.5,]
    con_id_z2 = QPN.add_constraint!(qp_net, cons_z2, lb_z2, ub_z2)

    cost_x = y-x
    level = 1
    QPN.add_qp!(qp_net, level, cost_x, [con_id_x,], x)

    cost_y = y
    level = 2
    QPN.add_qp!(qp_net, level, cost_y, [con_id_y], y)

    cost_z1 = (z[1]-x)^2
    level = 3
    QPN.add_qp!(qp_net, level, cost_z1, [con_id_z1,], z[1])
    
    cost_z2 = (z[2]-x)^2
    level = 3
    QPN.add_qp!(qp_net, level, cost_z2, [con_id_z2,], z[2])
    
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end
