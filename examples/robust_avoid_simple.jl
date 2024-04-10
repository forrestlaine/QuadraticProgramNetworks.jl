function setup(::Val{:robust_avoid_simple};
                 num_obj=2,
                 num_poly_faces=5,
                 exploration_vertices=10,
                 max_ego_delta=15.0,
                 max_obj_delta=1.0,
                 num_projections=5,
                 rng = MersenneTwister(1),
                 max_accel = 10.0,
                 kwargs...)

    # Ego poly: x1ₜ[1:2] + { z : Aₑz + bₑ ≥ 0 }
    # Obj poly: x2ₜ[1:2] + { z : Aₒz + bₒ ≥ 0 }
    
    Q = [0.0 0 ;
         0 0.001]
    q = [-1.0, 0]
    R = zeros(2,2)

    Ae_angles = collect(0:2*pi/num_poly_faces:1.99*pi) + 0.15*randn(rng,num_poly_faces) .+ pi*rand(rng)
    Ae = vcat(([cos(ai), sin(ai)]' for ai in Ae_angles)...)
    be = 0.2 .+ 0.8*rand(rng)*ones(num_poly_faces)

    Aos = map(1:num_obj) do i
        Ao_angles = collect(0:2*pi/num_poly_faces:1.99*pi) + 0.15*randn(rng,num_poly_faces) .+ pi*rand(rng)
        Ao_i = vcat(([cos(ai), sin(ai)]' for ai in Ao_angles)...)
    end
    bos = [0.2 .+ 0.8*rand(rng)*ones(num_poly_faces) for _ in 1:num_obj]

    ue = QPN.variables(:ue, 1:2)  
    uo = QPN.variables(:uo, 1:2, 1:num_obj)
    xe = QPN.variables(:xe, 1:2)  
    xo = QPN.variables(:xo, 1:2, 1:num_obj)

    s = QPN.variables(:s, 1:2, 1:num_obj)
    ϵ = QPN.variables(:ϵ, 1:num_obj)
 
    qp_net = QPNet(xe,xo, ue, uo, s, ϵ)
    qp_net.problem_data[:Ae] = Ae
    qp_net.problem_data[:be] = be
    qp_net.problem_data[:Ao] = Aos
    qp_net.problem_data[:bo] = bos
   
    s_players = Dict()

    for i = 1:num_obj
        cost = ϵ[i]
        cons = [Ae*(s[:,i]-(xe+ue)) + be + ones(num_poly_faces)*ϵ[i]; 
                Aos[i]*(s[:,i]-(xo[:,i]+uo[:,i])) + bos[i] + ones(num_poly_faces)*ϵ[i]]
        lb = fill(0.0, length(cons))
        ub = fill(Inf, length(cons))
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        vars = [s[:,i]; ϵ[i]]
        s_players[i] = QPN.add_qp!(qp_net, cost, [con_id,], vars...)
    end
 
    a_players = Dict()

    for i = 1:num_obj
        cons = uo[:,i]
        lb = fill(-max_obj_delta, 2)
        ub = fill(+max_obj_delta, 2)
        ad_con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        cost = ϵ[i]
        a_players[i] = QPN.add_qp!(qp_net, cost, [ad_con_id], uo[:,i])
    end

    cons = [ue; ϵ]
    lb = [fill(-max_ego_delta, 2); zeros(num_obj)]
    ub = [fill(+max_ego_delta, 2); fill(Inf, num_obj)]
    ego_con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
    xef = xe+ue
    cost = 0.5 * xef'*Q*xef + xef'*q + 0.5*ue'*R*ue
    ego_player = QPN.add_qp!(qp_net, cost, [ego_con_id], ue)

    #################
    # Add edges
    edge_list = [[(ego_player, a_players[i]) for i in 1:num_obj]; [(a_players[i], s_players[i]) for i in 1:num_obj]]

    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=exploration_vertices, num_projections=num_projections, debug_visualize=true, kwargs...)

    x0e = [-5.0, 0]
    x0os = map(1:num_obj) do i
        [3*(i-1), -1.0]
    end
    init = [x0e; vcat(x0os...); zeros(2*(1+num_obj)); zeros(3*num_obj)]
    qp_net.default_initialization .= init
    qp_net.visualization_function=(x->viz_solution_robust_avoid_simple(Ae,be,Aos,bos,x))

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

function viz_solution_robust_avoid_simple(Ae,be,Ao,bo,θ)
    f = Figure()
    ax = Axis(f[1,1], aspect = DataAspect(), limits = (-7, 15, -4, 18))
    #ax = Axis(f[1,1], aspect = DataAspect())
    num_obj = length(Ao)

    Xe = θ[1:2]
    Xo = [θ[2+(i-1)*2+1:2+i*2] for i in 1:num_obj]
    Ue = θ[(1+num_obj)*2+1:(1+num_obj)*2 + 2]
    Uo = [θ[(2+num_obj)*2+(i-1)*2+1:(2+num_obj)*2+i*2] for i in 1:num_obj]
    ϵ = θ[end-num_obj+1:end]

    Ve = get_verts(Ae,be,Xe)
    xe = [[v[1] for v in Ve]; Ve[1][1]]
    ye = [[v[2] for v in Ve]; Ve[1][2]]
    lines!(ax, xe, ye, color=:blue)
    Ve = get_verts(Ae,be,Xe+Ue)
    xe = [[v[1] for v in Ve]; Ve[1][1]]
    ye = [[v[2] for v in Ve]; Ve[1][2]]
    lines!(ax, xe, ye, color=:blue, linewidth=3)

    colors = [:red, :orange, :purple, :yellow]
    for i in 1:num_obj
        Vo = get_verts(Ao[i],bo[i],Xo[i])
        xo = [[v[1] for v in Vo]; Vo[1][1]]
        yo = [[v[2] for v in Vo]; Vo[1][2]]

        lines!(ax, xo, yo, color=colors[i])
        Vo = get_verts(Ao[i],bo[i],Xo[i]+Uo[i])
        xo = [[v[1] for v in Vo]; Vo[1][1]]
        yo = [[v[2] for v in Vo]; Vo[1][2]]
        lines!(ax, xo, yo, color=colors[i], linewidth=3)
    
        Ve = get_verts(Ae,be .+ ϵ[i],Xe+Ue)
        xe = [[v[1] for v in Ve]; Ve[1][1]]
        ye = [[v[2] for v in Ve]; Ve[1][2]]
        lines!(ax, xe, ye, color=colors[i], linestyle=:dash, linewidth=2)
    
        Vo = get_verts(Ao[i],bo[i] .+ ϵ[i], Xo[i]+Uo[i])
        xo = [[v[1] for v in Vo]; Vo[1][1]]
        yo = [[v[2] for v in Vo]; Vo[1][2]]
        lines!(ax, xo, yo, color=colors[i], linestyle=:dash, linewidth=2)
    end
    name = string(floor(Int, time()*1e3))[end-4:end] * ".png"
    display(f)
    f
end
