function simple_dyn(x, u)
    x + [x[3:4]; u]
end

function setup(::Val{:robust_avoid_simple}; T=1,
                 num_obj=1,
                 num_obj_faces=4,
                 obstacle_spacing = 1.0,
                 exploration_vertices=0,
                 num_projections=5,
                 lane_heading = 0.0,
                 initial_speed=3.0,
                 lane_width = 10.0,
                 initial_box_length = 6.0,
                 lane_dist_incentive = 10.0,
                 max_accel = 10.0,
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

    Ae = [1.0 1.1;
          -1  1.1;
          -1.2 -1;
          1 -0.85]
    be = ones(4)

    bo = 0.3*ones(3)
    Ao = [1.0 0.1;
          -1 0.1;
          0.2 1;
          -0.1 -1]
    bo = [0.2,0.2,0.6,0.6]

    x1 = QPN.variables(:x1, 1:4, 1:T)  
    x2 = QPN.variables(:x2, 1:4, 1:T)

    u1 = QPN.variables(:u1, 1:2, 1:T)
    u2 = QPN.variables(:u2, 1:2, 1:T)
    x̄1 = QPN.variables(:x̄1, 1:4)
    x̄2 = QPN.variables(:x̄2, 1:4)

    s = QPN.variables(:s, 1:2, 1:T)
    ϵ = QPN.variables(:ϵ, 1:T)
 
    qp_net = QPNet(x̄1,x̄2,x1,x2,u1,u2,s,ϵ)
    qp_net.problem_data[:Ae] = Ae
    qp_net.problem_data[:be] = be
    qp_net.problem_data[:Ao] = Ao
    qp_net.problem_data[:bo] = bo
    qp_net.problem_data[:T] = T
   
    #####################################################################
    
    # Add players responsible for identifying least-violated obstacle halfspace
    # constraint (min s s.t. s ≥ (aᵢ'x - bᵢ)) 
    # [recall obstacle avoidance ⇿ ∃ i s.t. aᵢ'x > bᵢ, i.e. only one halfspace
    # constraint needs to be satisfied to "avoid" the obstacle]
    
    s_players = Dict()

    for t = 1:T
        cost = ϵ[t]
        cons = [Ae*(s[:,t]-x1[1:2,t]) + be + ones(4)*ϵ[t]; Ao*(s[:,t]-x2[1:2,t]) + bo + ones(length(bo))*ϵ[t]]
        lb = fill(0.0, length(cons))
        ub = fill(Inf, length(cons))
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        vars = [s[:,t]; ϵ[t]]
        s_players[t] = QPN.add_qp!(qp_net, cost, [con_id,], vars...)
    end
  
    # Adversary
    dynamic_cons = []
    control_cons = []
    for t = 1:T
        if t==1
            append!(dynamic_cons, x2[:, t] - simple_dyn(x̄2, u2[:, t]))
        else
            append!(dynamic_cons, x2[:, t] - simple_dyn(x2[:, t-1], u2[:, t]))
        end
        append!(control_cons, u2[:,t])
    end
    lb_dyn = fill(0.0, length(dynamic_cons))
    ub_dyn = fill(0.0, length(dynamic_cons))
    lb_control = fill(-0.00, length(control_cons))
    ub_control = fill(+0.00, length(control_cons))
    ad_con_id = QPN.add_constraint!(qp_net, [dynamic_cons; control_cons], [lb_dyn; lb_control], [ub_dyn; ub_control])
    cost = sum(ϵ[end])
    #cost = sum(x2[1:2,t]'*x2[1:2,t] for t in 1:T)
    ad_player = QPN.add_qp!(qp_net, cost, [ad_con_id], x2, u2)

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
    lb_control = fill(-1.0, length(control_cons))
    ub_control = fill(+1.0, length(control_cons))
    ego_con_id = QPN.add_constraint!(qp_net, [dynamic_cons; control_cons], [lb_dyn; lb_control], [ub_dyn; ub_control])
    avoid_con_id = QPN.add_constraint!(qp_net, ϵ, zeros(T), fill(Inf, T))
    cost = sum(0.5*x1[:,t]'*Q*x1[:,t] + x1[:,t]'*q + 0.5*u1[:,t]'R*u1[:,t] for t in 1:T)
    ego_player = QPN.add_qp!(qp_net, cost, [ego_con_id, avoid_con_id], x1, u1)
    #ego_player = QPN.add_qp!(qp_net, cost, [ego_con_id, ], x1, u1)

    #################
    # Add edges
    edge_list = [(ego_player, ad_player); [(ad_player, s_player) for s_player in values(s_players)]]
    
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=exploration_vertices, num_projections=num_projections, visualization_function=(x->viz_solution(Ae,be,Ao,bo,x,T)), debug_visualize=true, kwargs...)
    
    qp_net
end

function get_verts(A,b,x)
    V = []
    for i in 1:length(b)
        for j in i+1:length(b)
            M = A[[i;j],:]
            m = b[[i;j]] - M*x
            try
                v = -M\m
                if all(A*(v-x) + b .≥ -1e-4)
                    push!(V,v)
                end
            catch err
            end
        end
    end
    c = sum(V) / length(V)
    angles = map(1:length(V)) do i
        d = V[i] - c
        atan(d[2],d[1])
    end
    I = sortperm(angles) |> reverse
    V[I]
end

function initialize(qpn, x0)
    θ = zeros(length(qpn.variables))
    θ[1:8] = x0

    x1prev = x0[1:4]
    x2prev = x0[5:8]
    T = qpn.problem_data[:T] 
    for t in 1:T
        u1 = [0.5,0]
        u2 = [0.0,0]
        x1prev = simple_dyn(x1prev, u1)
        θ[8+(t-1)*4+1:8+4*t] = x1prev
        x2prev = simple_dyn(x2prev, u2)
        θ[8+T*4+(t-1)*4+1:8+T*4+4*t] = x2prev
    end
    θ
end

function viz_solution(Ae,be,Ao,bo,θ,T)
    f = Figure()
    ax = Axis(f[1,1], aspect = DataAspect())


    Xe = θ[1:2]
    Xo = θ[5:6]

    Ve = get_verts(Ae,be,Xe)
    Vo = get_verts(Ao,bo,Xo)

    xe = [[v[1] for v in Ve]; Ve[1][1]]
    ye = [[v[2] for v in Ve]; Ve[1][2]]
    xo = [[v[1] for v in Vo]; Vo[1][1]]
    yo = [[v[2] for v in Vo]; Vo[1][2]]

    lines!(ax, xe, ye, color=:blue)
    lines!(ax, xo, yo, color=:red)

    for t in T:T
        Xe = θ[8+(t-1)*4+1:8+(t-1)*4+4]
        Xo = θ[8+T*4+(t-1)*4+1:8+T*4+(t-1)*4+4]
        ϵ = θ[end-T+t]

        Ve = get_verts(Ae,be,Xe[1:2])
        Vo = get_verts(Ao,bo,Xo[1:2])
        
        Ve_i = get_verts(Ae,be.+ϵ,Xe[1:2])
        Vo_i = get_verts(Ao,bo.+ϵ,Xo[1:2])

        xe = [[v[1] for v in Ve]; Ve[1][1]]
        ye = [[v[2] for v in Ve]; Ve[1][2]]
        xo = [[v[1] for v in Vo]; Vo[1][1]]
        yo = [[v[2] for v in Vo]; Vo[1][2]]
        lines!(ax, xe, ye, color=:blue, linewidth=t)
        lines!(ax, xo, yo, color=:red, linewidth=t)
        
        xe_i = [[v[1] for v in Ve_i]; Ve_i[1][1]]
        ye_i = [[v[2] for v in Ve_i]; Ve_i[1][2]]
        xo_i = [[v[1] for v in Vo_i]; Vo_i[1][1]]
        yo_i = [[v[2] for v in Vo_i]; Vo_i[1][2]]
        lines!(ax, xe_i, ye_i, color=:blue, linewidth=t, linestyle=:dash)
        lines!(ax, xo_i, yo_i, color=:red, linewidth=t, linestyle=:dash)
    end
    display(f)
end
