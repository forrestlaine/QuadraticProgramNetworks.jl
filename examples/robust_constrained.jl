using Symbolics
using StaticArrays
using Infiltrator
using QPN

struct PolyObject{D, T}
    # Assumed clockwise ordering
    verts::Vector{SVector{D, T}}
end

struct HalfSpace{D, T}
    a::SVector{D, T} # normal
    b::T             # offset
end

function halfspaces(poly_obj::PolyObject{D, T}) where {D, T}
    N = length(poly_obj.verts)
    map(1:N) do i
        ii = mod1(i+1, N)
        d = poly_obj.verts[ii] - poly_obj.verts[i]
        a = SVector(-d[2], d[1])
        b = a'*poly_obj.verts[i]
        HalfSpace(a, b)
    end
end

function dyn(x, u; Δ = 0.1)
    x + Δ * [x[3:4]+0.5*Δ*u[1:2]; 
             u[1:2]]
end

function setup(; T=2,
                 num_obj=1,
                 num_obj_faces=5,
                 obstacle_spacing = 4.0,
                 lane_heading = 0.0,
                 max_initial_speed=2.0,
                 lane_width = 10.0,
                 initial_box_length = 6.0,
                 lane_dist_incentive = 10.0,
                 max_accel = 1.0,
                 kwargs...)
        

    lane_vec = [cos(lane_heading), sin(lane_heading)]
    
    x = Symbolics.variables(:x, 1:T, 1:4)
    u = Symbolics.variables(:u, 1:T, 1:2)
    x̄ = Symbolics.variables(:x̄, 1:4)
    s = Symbolics.variables(:s, 1:T, 1:num_obj)
    o = Symbolics.variables(:o, 1:num_obj, 1:2)
    c₋ = Symbolics.variable(:c₋)
    
    qp_net = QPNet(x,u,x̄,s,o,c₋)
    
    objs = map(1:num_obj) do i
        verts = map(1:num_obj_faces) do j
            θ = j*2π / num_obj_faces    
            SVector(o[i,1]+cos(θ), o[i,2]+sin(θ))
        end
        PolyObject(verts)
    end
    obstacle_distances_along = collect(1:num_obj) * obstacle_spacing .+ initial_box_length / 2
    obstacle_offsets = [-1^i for i in 1:num_obj] .* lane_width/4.0
    
    obj_halfspaces = [halfspaces(obj) for obj in objs]
    right_laneline_normal = SVector(-sin(lane_heading), cos(lane_heading))
    lane_halfspaces = [HalfSpace(right_laneline_normal, -lane_width/2), 
                       HalfSpace(-right_laneline_normal, -lane_width/2)]
   
    #####################################################################
    
    # Add players responsible for identifying least-violated obstacle halfspace
    # constraint (min s s.t. s ≥ (aᵢ'x - bᵢ)) 
    # [recall obstacle avoidance ⇿ ∃ i s.t. aᵢ'x > bᵢ]
    for t = 1:T
        for i = 1:num_obj
            cost = s[t,i] # min s[t,i]
            cons = map(1:num_obj_faces) do j
                hs = obj_halfspaces[i][j]
                s[t,i] - (hs.a'*x[t,1:2] - hs.b)
            end
            lb = fill(0.0, num_obj_faces)
            ub = fill(Inf, num_obj_faces)
    
            con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
            
            level = 3
            player_id = QPN.add_qp!(qp_net, level, cost, [con_id,], s[t,i])
        end
    end
    
    #####################################################################
   
    # Add player responsible for finding most-violated constraint
    # (max c s.t. c <= all other constraints)
    min_cons = []
    for t = 1:T
        for lane_hs in lane_halfspaces
            push!(min_cons, (lane_hs.a'*x[t,1:2] - lane_hs.b) - c₋)
        end
        for i in 1:num_obj
            push!(min_cons, s[t,i] - c₋)
        end
    end
    lb = fill(0.0, length(min_cons))
    ub = fill(Inf, length(min_cons))
    con_id = QPN.add_constraint!(qp_net, min_cons, lb, ub)
    level = 3
    cost = -c₋
    player_id = QPN.add_qp!(qp_net, level, cost, [con_id,], c₋)

    #####################################################################

    # Add player responsible for choosing initial state and obstacle centers
    # to draw trajectory to boundary of infeasibility
   
    # Setup Dynamic Equation constraints
    dynamic_cons = []
    for t = 1:T
        if t==1
            append!(dynamic_cons, x[t,:] - dyn(x̄, u[t,:]))
        else
            append!(dynamic_cons, x[t,:] - dyn(x[t-1,:], u[t,:]))
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
    lb = -[initial_box_length/2, lane_width/2, max_initial_speed, max_initial_speed]
    ub = [initial_box_length/2, lane_width/2, max_initial_speed, max_initial_speed]
    init_con_id = QPN.add_constraint!(qp_net, initial_state_cons, lb, ub)

    # Setup Obstacle Constraints
    obstacle_cons = []
    lb = []
    ub = []
    for i in 1:num_obj
        append!(obstacle_cons, R\o[i,:])
        append!(lb, [obstacle_distances_along[i], obstacle_offsets[i]-lane_width/4])
        append!(ub, [obstacle_distances_along[i], obstacle_offsets[i]+lane_width/4])
    end
    obstacle_con_id = QPN.add_constraint!(qp_net, obstacle_cons, lb, ub)

    level = 2
    cost = c₋^2
    player_id = QPN.add_qp!(qp_net, level, cost, [dyn_con_id, init_con_id, obstacle_con_id], x̄, x, o)
    
     
    #####################################################################

    # Add player responsible for choosing control variables
    # to avoid worst-case obstacles, initial condition

    cons = [c₋, ]
    for t = 1:T
        append!(cons, u[t,:])
    end
    lb = [0; fill(-max_accel, 2*T)]
    ub = [Inf; fill(max_accel, 2*T)]
    con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
    
    primary_cost = sum(0.5*u[t,:]'*u[t,:] - lane_dist_incentive * x[t,1:2]'*lane_vec for t = 1:T)
    level = 1
    player_id = QPN.add_qp!(qp_net, level, primary_cost, [con_id,], u)

    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end
