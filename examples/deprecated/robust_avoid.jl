function setup(::Val{:robust_avoid}; T=3,
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
    c = QPN.variable(:c)
    
    qp_net = QPNet(o,x̄,x,u,h,s,c)
 
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
   
    #### Adversary
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
    

    # Setup Initial State constraints
    initial_state_cons = []
    R = [lane_vec right_laneline_normal]
    append!(initial_state_cons, R\x̄[1:2])
    append!(initial_state_cons, x̄[3:4])
    lb = [0, 0.0, initial_speed, 0.0]
    ub = [0, 0.0, initial_speed, 0.0]
    init_con_id = QPN.add_constraint!(qp_net, initial_state_cons, lb, ub)

    obstacle_cons = vcat(([o[1,i], o[2,i]] for i in 1:num_obj)...)
    lb = repeat([1.5, -1], num_obj)
    ub = repeat([2.5, 1], num_obj)
    o_con_id = QPN.add_constraint!(qp_net, obstacle_cons, lb, ub)
    
    cost = c
    a_player = QPN.add_qp!(qp_net, cost, [dyn_con_id, init_con_id, o_con_id], x̄, o, x)


    #### C player
    cs_cons = []
    for t in 1:T
        for i in 1:num_obj
            push!(cs_cons, c-s[i,t])
        end
    end
    lb = fill(-Inf, length(cs_cons))
    ub = fill(0.0, length(cs_cons))
    cs_con_id = QPN.add_constraint!(qp_net, cs_cons, lb, ub)

    cost = -c
    c_player = QPN.add_qp!(qp_net, cost, [cs_con_id], c)
    
    ##### protagonist

    control_cons = []
    for t = 1:T
        append!(control_cons, u[:, t])
    end
    lb = fill(-max_accel, 2*T)
    ub = fill(max_accel, 2*T)
    control_con_id = QPN.add_constraint!(qp_net, control_cons, lb, ub)

    cons = [s[i,t] for i in 1:num_obj for t in 1:T]
    lb = fill(0.0, length(cons))
    ub = fill(Inf, length(cons))
    s_con_id = QPN.add_constraint!(qp_net, cons, lb, ub)

    #cost = sum((u[1,t]-15)^2+u[2,t]^2 for t = 1:T)
    cost = -10*c+sum(-10*x[1,t]+1*x[2,t]^2 for t = 1:T)
    u_player = QPN.add_qp!(qp_net, cost, [control_con_id], u)

    #################
    # Add edges
    edge_list = [(u_player, a_player); (a_player, c_player); [(a_player, s_player) for s_player in values(s_players)]]
    
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=exploration_vertices, num_projections=num_projections, kwargs...)
    
    qp_net
end
