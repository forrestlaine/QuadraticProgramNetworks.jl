"""
variables := w1 w2 x y
f1: ½ || [x; y] - [w1; w2] ||²
f2: ½ (y-x)²
"""
function setup(::Val{:simple_bilevel}; kwargs...)
    
    w = Symbolics.variables(:w, 1:2) 
    x = Symbolics.variable(:x)
    y = Symbolics.variable(:y)
    
    qp_net = QPNet(w,x,y)
  
    lb = [0.0,]
    ub = [Inf,]
    cons = [y,]
    con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
           
    cost = (y-x)^2
    qp_id1 = QPN.add_qp!(qp_net, cost, [con_id,], y)

    cost = ([x;y] - w)'*([x; y] - w)
    qp_id2 = QPN.add_qp!(qp_net, cost, [], x)
   
    edge_list = [(qp_id2, qp_id1)]

    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    qp_net
end

