function simple_pointmass_dyn(x, u, dt)
    x + dt * [x[3:4] + dt/2 * u; u]
end

function setup(::Val{:shepherd_sheep}; 
               T=5,
               dt = 1.0,
               shepherd_start=[-4.0, -5, 0, 0],
               sheep_start=[5.0, -4, 0, 0],
               u_max = 1.0,
               kwargs...)

    X_shepherd = QPN.variables(:X_shepherd, 1:4, 1:T)
    X_sheep    = QPN.variables(:X_sheep   , 1:4, 1:T)
    U_shepherd = QPN.variables(:U_shepherd, 1:2, 1:T)
    U_sheep    = QPN.variables(:U_sheep   , 1:2, 1:T)

    qp_net = QPNet(X_shepherd, X_sheep, U_shepherd, U_sheep)

    running_shepherd_costs = []
    running_sheep_costs = []
 
    shepherd_players = Dict()
    sheep_players = Dict()
    shepherd_cons = []

    for t = T:-1:1
        # shepherd
        x1 = X_shepherd[:,t]
        u1 = U_shepherd[:,t]
        x2 = X_sheep[:,t]
        u2 = U_sheep[:,t]

        Xprev = t > 1 ? X_shepherd[:, t-1] : shepherd_start
        dyn_cons = x1 - simple_pointmass_dyn(Xprev, u1, dt)
        #cons = [u1; dyn_cons]
        #l = [fill(-u_max, 2); zeros(4)]
        #u = [fill(+u_max, 2); zeros(4)]
        cons = dyn_cons
        l = zeros(4)
        u = zeros(4)
        con_id = QPN.add_constraint!(qp_net, cons, l, u)
        push!(shepherd_cons, con_id)
        
        cost = x2'*x2 + 0.1*u1'*u1
        push!(running_shepherd_costs, cost)
        shepherd_players[t] = QPN.add_qp!(qp_net, sum(running_shepherd_costs), [con_id,], x1, u1)
        
        Xprev = t > 1 ? X_sheep[:, t-1] : sheep_start
        dyn_cons = x2 - simple_pointmass_dyn(Xprev, u2, dt)
        #cons = [u2; dyn_cons]
        #l = [fill(-u_max, 2); zeros(4)]
        #u = [fill(+u_max, 2); zeros(4)]
        cons = dyn_cons
        l = zeros(4)
        u = zeros(4)
        con_id = QPN.add_constraint!(qp_net, cons, l, u)
        
        cost = (x2-x1)'*(x2-x1) + 0.1*u2'*u2
        push!(running_sheep_costs, cost)
        sheep_players[t] = QPN.add_qp!(qp_net, sum(running_sheep_costs), [con_id,], x2, u2)

    end
    #shepherd_player = QPN.add_qp!(qp_net, sum(running_shepherd_costs), shepherd_cons, X_shepherd, U_shepherd)

    edge_list = []
    for t in 1:T-1
        #push!(edge_list, (shepherd_players[t], shepherd_players[t+1]))
        push!(edge_list, (shepherd_players[t], sheep_players[t]))
        push!(edge_list, (sheep_players[t], shepherd_players[t+1]))
        #push!(edge_list, (sheep_players[t], sheep_players[t+1]))
    end
    push!(edge_list, (shepherd_players[T], sheep_players[T]))
    
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=10, check_convexity=true, debug_visualize=true, kwargs...)
    qp_net.default_initialization = [repeat(shepherd_start, T); repeat(sheep_start, T); zeros(2*T*2)]
    qp_net.visualization_function=(x->visualize_shepherd_sheep(shepherd_start, sheep_start, T, x))

    qp_net
end

function visualize_shepherd_sheep(x1_0, x2_0, T, θ)
    f = Figure()
    ax = Axis(f[1,1], aspect = DataAspect(), limits = (-10, 10, -10, 10))

    scatter!(ax, x1_0[1], x1_0[2], color=:blue, markersize=10)
    scatter!(ax, x2_0[1], x2_0[2], color=:red, markersize=10)

    for t in 1:T
        x1 = θ[(t-1)*4+1:(t-1)*4+2]
        scatter!(ax, x1[1], x1[2], color=:blue, markersize=10)
        x2 = θ[T*4+(t-1)*4+1:T*4+(t-1)*4+2]
        scatter!(ax, x2[1], x2[2], color=:red, markersize=10)
    end
    display(f)
    f
end
