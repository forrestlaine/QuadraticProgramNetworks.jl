using Symbolics
using QPN


"""
variables := w1 w2 x y
f1: ½ || [x; y] - [w1; w2] ||²
f2: ½ (y-x)²
"""
function setup(; debug=false,
                 shared_variable_mode=QPN.SHARED_DUAL,
                 high_dimension=false,
                 gen_solution_map=true,
              )
    
    w = Symbolics.variables(:w, 1:2) 
    x = Symbolics.variable(:x)
    y = Symbolics.variable(:y)
    
    qp_net = QPNet(w,x,y)

  
    lb = [0.0,]
    ub = [Inf,]
    cons = [y,]
    con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
           
    cost = (y-x)^2
    level = 2
    QPN.add_qp!(qp_net, level, cost, [con_id,], y)

    cost = ([x;y] - w)'*([x; y] - w)
    level = 1
    QPN.add_qp!(qp_net, level, cost, [], x)
    
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; debug, shared_variable_mode, high_dimension, gen_solution_map)

    qp_net
end

