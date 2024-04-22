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
    con_id = QPNets.add_constraint!(qp_net, cons, lb, ub)
           
    cost = (y-x)^2
    qp_id1 = QPNets.add_qp!(qp_net, cost, [con_id,], y)

    cost = ([x;y] - w)'*([x; y] - w)
    qp_id2 = QPNets.add_qp!(qp_net, cost, [], x)
   
    edge_list = [(qp_id2, qp_id1)]

    QPNets.add_edges!(qp_net, edge_list)
    QPNets.assign_constraint_groups!(qp_net)
    QPNets.set_options!(qp_net; debug_visualize=false, kwargs...)
    qp_net.default_initialization .= zeros(4)
    
    qp_net.visualization_function=visualize_simple_bilevel

    qp_net
end


function visualize_simple_bilevel(θ)
    f = Figure()
    ax = Axis(f[1,1], aspect = DataAspect(), limits = (-4, 4, -3, 5))

    w = θ[1:2]
    x = θ[3]
    y = θ[4]

    lines!(ax, [-5,0.0], [0.0, 0], color=:black, linewidth=4)
    lines!(ax, [0.0, 5.0], [0.0, 5.0], color=:black, linewidth=4)
    
    scatter!(ax, w[1],w[2], color=:green, markersize=25)
    scatter!(ax, x,y, color=:blue, markersize=25)
    name = string(floor(Int, time()*1e3))[end-4:end] * ".png"
    save(name, f)
    display(f)
end



