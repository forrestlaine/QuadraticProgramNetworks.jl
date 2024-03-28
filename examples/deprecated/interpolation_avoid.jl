function setup(::Val{:interpolation_avoid}; 
               T=1,
               kwargs...)

    # Ego poly: x1ₜ[1:2] + { z : Aₑz + bₑ ≥ 0 }
    # Obj poly: x2ₜ[1:2] + { z : Aₒz + bₒ ≥ 0 }
    
    Q = [0.0 0 0 0 ;
         0 0 0 0;
         0 0 0 0;
         0 0 0 0]
    R = [0. 0;
         0 0.]
    q = [-1.0, 0, 0, 0]

    Aₑ = [1.0 1;
          -1  1;
          -1 -1;
          1 -1]
    bₑ = ones(4)

    a2 = [sqrt(3), -2]
    a2 = a2 / norm(a2)
    a3 = [-sqrt(3), -2]
    a3 = a3 / norm(a3)

    Aₒ = [0 1.0;
          a2'
          a3']
    bₒ = 0.3*ones(3)

    x1 = QPN.variables(:x1, 1:4, 1:T)  
    x2 = QPN.variables(:x2, 1:4, 1:T)

    u = QPN.variables(:u1, 1:2, 1:T)
    p = QPN.variables(:p, 1:T)
    x̄ = QPN.variables(:x̄, 1:4)

    s = QPN.variables(:s, 1:2, 1:T)
    ϵ = QPN.variables(:ϵ, 1:T)
 
    qp_net = QPNet(x̄,x1,x2,u,p,s,ϵ)
   
    s_players = Dict()
    a_players = Dict()

    for t = 1:T
        cost = ϵ[t]
        cons = [Aₑ*(s[:,t]-x2[1:2,t]) + bₑ + ones(4)*ϵ[t]; Aₒ*(s[:,t]) + bₒ + ones(3)*ϵ[t]]
        lb = fill(0.0, length(cons))
        ub = fill(Inf, length(cons))
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        vars = [s[:,t]; ϵ[t]]
        s_players[t] = QPN.add_qp!(qp_net, cost, [con_id,], vars...)

        cost = ϵ[t]
        if t == 1
            cons = [x2[:,t] - p[t]*x̄[:] - (1.0-p[t])*x1[:,t]; p[t]]
        else
            cons = [x2[:,t] - p[t]*x1[:,t-1] - (1.0-p[t])*x1[:,t]; p[t]]
        end
        lb = fill(0.0, length(cons))
        ub = [fill(0.0, length(cons)-1); 1.0]
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        vars = [x2[:,t]; p[t]]
        a_players[t] = QPN.add_qp!(qp_net, cost, [con_id], vars...)
    end

    # Ego
    dynamic_cons = []
    control_cons = []
    for t = 1:T
        if t==1
            append!(dynamic_cons, x1[:, t] - simple_dyn(x̄1, u1[:, t]))
        else
            append!(dynamic_cons, x1[:, t] - simple_dyn(x1[:, t-1], u1[:, t]))
        end
        append!(control_cons, u1[:,t])
    end
    lb_dyn = fill(0.0, length(dynamic_cons))
    ub_dyn = fill(0.0, length(dynamic_cons))
    lb_control = fill(-5.0, length(control_cons))
    ub_control = fill(+5.0, length(control_cons))
    ego_con_id = QPN.add_constraint!(qp_net, [dynamic_cons; control_cons; ϵ], [lb_dyn; lb_control; zeros(T)], [ub_dyn; ub_control; fill(Inf, T)])
    cost = sum(0.5*x1[:,t]'*Q*x1[:,t] + x1[:,t]'*q + 0.5*u1[:,t]'R*u1[:,t] for t in 1:T)
    ego_player = QPN.add_qp!(qp_net, cost, [ego_con_id], x1, u1)

    edge_list = [[(ego_player, a_player) for a_player in values(a_players)]; [(a_players[t], s_players[t]) for t in 1:T]]
    
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=exploration_vertices, num_projections=num_projections, kwargs...)
    
    qp_net, Aₑ, bₑ, Aₒ, bₒ
end

function initialize_interpolation(qpn, x0, T)
    θ = zeros(length(qpn.variables))
    θ[1:8] = x0

    x1prev = x0[1:4]
    x2prev = x0[5:8]
    
    for t in 1:T
        u = zeros(2)
        x1prev = simple_dyn(x1prev, u)
        θ[8+(t-1)*4+1:8+4*t] = x1prev
        x2prev = simple_dyn(x2prev, u)
        θ[8+T*4+(t-1)*4+1:8+T*4+4*t] = x2prev
    end
    θ
end

function viz_interpolation_solution(Ae,be,Ao,bo,θ,T)
    f = Figure(; aspect=DataAspect())
    ax = Axis(f[1,1])


    Xe = θ[1:2]
    Xo = θ[5:6]

    Ve = get_verts(Ae,be,Xe)
    Vo = get_verts(Ao,bo,Xo)

    xe = [[v[1] for v in Ve]; Ve[1][1]]
    ye = [[v[2] for v in Ve]; Ve[1][2]]
    xo = [[v[1] for v in Vo]; Vo[1][1]]
    yo = [[v[2] for v in Vo]; Vo[1][2]]

    lines!(ax, xe, ye, color=:blue, linestyle=:dash)
    lines!(ax, xo, yo, color=:red, linestyle=:dash)

    for t in 1:T
        Xe = θ[8+(t-1)*4+1:8+(t-1)*4+2]
        Xo = θ[8+T*4+(t-1)*4+1:8+T*4+(t-1)*4+2]
        Ve = get_verts(Ae,be,Xe)
        Vo = get_verts(Ao,bo,Xo)

        xe = [[v[1] for v in Ve]; Ve[1][1]]
        ye = [[v[2] for v in Ve]; Ve[1][2]]
        xo = [[v[1] for v in Vo]; Vo[1][1]]
        yo = [[v[2] for v in Vo]; Vo[1][2]]
        lines!(ax, xe, ye, color=:blue, linewidth=2*t)
        lines!(ax, xo, yo, color=:red, linewidth=2*t)
    end
        

    display(f)
end
