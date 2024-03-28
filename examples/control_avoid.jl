function visualize(qpn, x; num_obj_faces=4, lane_width = 10.0, num_obj=1, T=3)
    f = Figure()
    f = Figure(resolution=(1440,1440))
    ax = f[1, 1] = Axis(f, aspect = DataAspect())
    xlims!(ax, -4.0, 12.0)
    ylims!(ax, -8.0, 8.0)

    lines!(ax, [-4.0, 12], [-lane_width/2, -lane_width/2], color=:black)
    lines!(ax, [-4.0, 12], [lane_width/2, lane_width/2], color=:black)

    vals = (; O = Dict(1=>x[1:2]), X0 = x[3:6], X = Dict(t=>x[6+(t-1)*4+1:6+4*t] for t in 1:T))

    p = Circle(Point(vals.X0[1:2]...), 0.1f0)
    scatter!(ax, p, color=:green)

    for i = 1:num_obj
        verts = map(1:num_obj_faces) do j
            θ = j*2π / num_obj_faces    
            SVector(vals.O[i][1]+cos(θ), vals.O[i][2]+sin(θ))
        end
        push!(verts, verts[1])
        lines!(ax, verts, color=:red)
    end

    for t = 1:T
        xt = vals.X[t]
        p = Circle(Point(xt[1:2]...), 0.1f0)
        scatter!(ax, p, color=:blue)
    end
    display(f) 
end

function setup(::Val{:control_avoid}; T=3,
                 num_obj=1,
                 num_obj_faces=4,
                 obstacle_spacing = 1.0,
                 exploration_vertices=10,
                 num_projections=5,
                 lane_heading = 0.0,
                 initial_speed=3.0,
                 lane_width = 10.0,
                 initial_box_length = 6.0,
                 lane_dist_incentive = 10.0,
                 max_accel = 10.0,
                 kwargs...)

    lane_vec = [cos(lane_heading), sin(lane_heading)]
    
    o = QPN.variables(:o, 1:2, 1:num_obj)
    x = QPN.variables(:x, 1:4, 1:T)
    u = QPN.variables(:u, 1:2, 1:T)
    x̄ = QPN.variables(:x̄, 1:4)
    s = QPN.variables(:s, 1:num_obj, 1:T)
    h = QPN.variables(:h, 1:num_obj_faces, 1:num_obj, 1:T)
    #c = QPN.variable(:c)
    
    qp_net = QPNet(o,x̄,x,u,h,s)
 
    objs = map(1:num_obj) do i
        verts = map(1:num_obj_faces) do j
            θ = j*2π / num_obj_faces    
            SVector(o[1, i]+cos(θ), o[2, i]+sin(θ))
        end
        PolyObject(verts)
    end
    obstacle_distances_along = collect(1:num_obj) * obstacle_spacing .+ initial_box_length / 2
    obstacle_offsets = [(-1)^i for i in 1:num_obj] .* lane_width/5.0
    
    obj_halfspaces = [halfspaces(obj) for obj in objs]
    right_laneline_normal = SVector(-sin(lane_heading), cos(lane_heading))
    lane_halfspaces = [HalfSpace(right_laneline_normal, -lane_width/2), 
                       HalfSpace(-right_laneline_normal, -lane_width/2)]
   
    #####################################################################
    
    # Add players responsible for identifying least-violated obstacle halfspace
    # constraint (min s s.t. s ≥ (aᵢ'x - bᵢ)) 
    # [recall obstacle avoidance ⇿ ∃ i s.t. aᵢ'x > bᵢ, i.e. only one halfspace
    # constraint needs to be satisfied to "avoid" the obstacle]
    
    s_players = Dict()

    for t = 1:T
        for i = 1:num_obj
            cost = s[i, t] # min s[t,i]
            cons = mapreduce(vcat, 1:num_obj_faces) do j
                hs = obj_halfspaces[i][j]
                [h[j, i, t] - (hs.a'*x[1:2,t] - hs.b),
                 s[i, t] - h[j, i, t]]
            end
            lb = fill(0.0, 2*num_obj_faces)
            ub = mapreduce(vcat, 1:num_obj_faces) do j
                [0.0, Inf]
            end
            con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
            
            vars = [s[i, t]; [h[j, i, t] for j in 1:num_obj_faces]]
            s_players[t,i] = QPN.add_qp!(qp_net, cost, [con_id,], vars...)
        end
    end
    
    #####################################################################
   
    # Add player responsible for finding most-violated constraint
    # (max c s.t. c <= all other constraints)
    #min_cons = []
    #for t = 1:T
    #    #for lane_hs in lane_halfspaces
    #    #    push!(min_cons, (lane_hs.a'*x[1:2,t] - lane_hs.b) - c)
    #    #end
    #    for i in 1:num_obj
    #        push!(min_cons, s[i,t] - c)
    #    end
    #end
    #lb = fill(0.0, length(min_cons))
    #ub = fill(Inf, length(min_cons))
    #con_id = QPN.add_constraint!(qp_net, min_cons, lb, ub)
    #cost = -c
    #c_player = QPN.add_qp!(qp_net, cost, [con_id,], c)

    #####################################################################

    # Add player responsible for driving 
    # to avoid obstacles
   
    dynamic_cons = []
    for t = 1:T
        if t==1
            append!(dynamic_cons, x[:, t] - dyn(x̄, u[:, t]))
        else
            append!(dynamic_cons, x[:, t] - dyn(x[:, t-1], u[:, t]))
        end
    end
    lb = fill(0.0, length(dynamic_cons))
    ub = fill(0.0, length(dynamic_cons))
    dyn_con_id = QPN.add_constraint!(qp_net, dynamic_cons, lb, ub)
    
    control_cons = []
    for t = 1:T
        append!(control_cons, u[:, t])
    end
    lb = fill(-max_accel, 2*T)
    ub = fill(max_accel, 2*T)
    control_con_id = QPN.add_constraint!(qp_net, control_cons, lb, ub)

    # Setup Initial State constraints
    initial_state_cons = []
    R = [lane_vec right_laneline_normal]
    append!(initial_state_cons, R\x̄[1:2])
    append!(initial_state_cons, x̄[3:4])
    lb = [0, 0.0, initial_speed, 0.0]
    ub = [0, 0.0, initial_speed, 0.0]
    init_con_id = QPN.add_constraint!(qp_net, initial_state_cons, lb, ub)

    cons = [s[i,t] for i in 1:num_obj for t in 1:T]
    lb = fill(0.0, length(cons))
    ub = fill(Inf, length(cons))
    s_con_id = QPN.add_constraint!(qp_net, cons, lb, ub)

    #cost = sum((u[1,t]-15)^2+u[2,t]^2 for t = 1:T)
    cost = sum(-10*x[1,t]+1*x[2,t]^2 for t = 1:T)
    u_player = QPN.add_qp!(qp_net, cost, [dyn_con_id, control_con_id, init_con_id, s_con_id], x̄, x, u)

    #################
    # Add edges

    #edge_list = [(u_player, c_player); [(c_player, s_player) for s_player in values(s_players)]]
    #edge_list = [(c_player, s_player) for s_player in values(s_players)]
    edge_list = [(u_player, s_player) for s_player in values(s_players)]
    
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=exploration_vertices, num_projections=num_projections, kwargs...)
    
    qp_net
end
