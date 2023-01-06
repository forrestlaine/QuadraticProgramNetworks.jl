using GLMakie
using GeometryBasics

function compute_vertices(p, normals, lengths)
    verts = Point[]
    for i in 1:length(normals)
        if i == 1
            A = [normals[end]'; normals[1]']
            b = [lengths[end]; lengths[1]] + A*p
        else
            A = [normals[i-1]'; normals[i]']
            b = [lengths[i-1]; lengths[i]] + A*p
        end
        vert = A\b
        push!(verts, Point(vert...))
    end
    verts
end
            
@testset "Bouncing Polygon" begin
    # Uncontrolled, no friction version
    ###################################
    
    # Parameters
    T = 15   # number of simulation steps
    Δ = 0.1  # simulation timestep
    p0 = [0.0; 3.5] # initial configuration
    v0 = [0; -1.5]  # initial velocity
    g = [0.0, 0.0]  # gravity vector
    Kₚ = 100.0 
    Kᵥ = 1.0
    M = 5.0 # centroid mass
    m = 0.5 # face mass

    num_faces = 4
    d = [[cos(θ), sin(θ)] for θ in collect(1:num_faces)*(2*pi)/num_faces]
    b̄ = 2.0*ones(Float64, num_faces)
    surface_normals = [[1.0, 3.0], [-1.0, 1.0], [1.0, -1.0], [-1.0, -3.0]]
    surface_nominals = [3.5, 0, -8.5, -20.0]
    num_surfaces = length(surface_nominals)
    
    # vars: [ p0 v0 p₁0 v₁0 p₂0 v₂0 ... pₙ0 vₙ0 | λ1 p1 v1 p₁1 v₁1 p₂1 v₂1 ... | λ2 p2 v2 ... | ... ]
   
    sim_state_dim = num_faces * 5 + 4
    total_dim = (num_faces+1) * 4 + T*sim_state_dim # add initial values into state for convenience

    qps = Dict{Int, QP}()
    sets =  Dict{Int, QPN.Poly}()
    level_1_progs = Set{Int}()
    level_2_progs = Set{Int}()
    dyn_sets = []
    
    x0 = [p0; v0; reduce(vcat, ([p0 + b̄[f]*d[f]; v0] for f in 1:num_faces))]
    x = [x0; reduce(vcat, ([zeros(num_faces); x0] for _ in 1:T))]
    
    qp_ind = 1
    set_ind = 1

    initial_offset = 4+num_faces*2

    for t ∈ 1:T
        A = zeros(4 * (num_faces+1), total_dim)
        l = zeros(4 * (num_faces+1))
        u = zeros(4 * (num_faces+1))

        offset = (num_faces+1)*4 + (t-1)*sim_state_dim + num_faces # start at | pt vt ... | part of sim state
        poffset = offset - sim_state_dim # previous state offset
        
        # p' - p - Δ v' = 0
        A[1, offset+1] = 1.0
        A[2, offset+2] = 1.0
        A[1, offset+3] = -Δ
        A[2, offset+4] = -Δ
        A[1, poffset+1] = -1.0
        A[2, poffset+2] = -1.0

        # fᵢ = Kₚ( dᵢ'(pᵢ-p)) + Kᵥ(dᵢ'(vᵢ-v)) - Kₚ*b̄ 
        #
        # v' - v - Δ( ∑fᵢ*dᵢ) / M = Δg 
        # v' - v - Δ/M * ∑ ( dᵢ* (Kₚ(dᵢ'(pᵢ-p)) + Kᵥ(dᵢ'(vᵢ-v)))) = Δg-Δ/M∑(dᵢKₚb̄)
        A[3, offset+3] = 1.0 
        A[3, poffset+3] = -1.0
        A[4, offset+4] = 1.0
        A[4, poffset+4] = -1.0
        l[3:4] = u[3:4] = Δ * g

        for f ∈ 1:num_faces
            A[3:4, poffset+1:poffset+2] += Kₚ*Δ/M * d[f] * d[f]'
            A[3:4, poffset+3:poffset+4] += Kᵥ*Δ/M * d[f] * d[f]'
            A[3:4, poffset+4+(f-1)*4+1:poffset+4+(f-1)*4+2] -= Kₚ*Δ/M * d[f]*d[f]'
            A[3:4, poffset+4+(f-1)*4+3:poffset+4+(f-1)*4+4] -= Kᵥ*Δ/M * d[f]*d[f]'

            l[3:4] -= Kₚ*Δ/M * d[f] * b̄[f]
            u[3:4] -= Kₚ*Δ/M * d[f] * b̄[f]

            A[4+(f-1)*4+1,offset+4+(f-1)*4+1] = 1.0
            A[4+(f-1)*4+1,poffset+4+(f-1)*4+1] = -1.0
            A[4+(f-1)*4+1,offset+4+(f-1)*4+3] = -Δ
            A[4+(f-1)*4+2,offset+4+(f-1)*4+2] = 1.0
            A[4+(f-1)*4+2,poffset+4+(f-1)*4+2] = -1.0
            A[4+(f-1)*4+2,offset+4+(f-1)*4+4] = -Δ

            A[4+(f-1)*4+3,offset+4+(f-1)*4+3] = 1.0
            A[4+(f-1)*4+3,poffset+4+(f-1)*4+3] = -1.0
            A[4+(f-1)*4+4,offset+4+(f-1)*4+4] = 1.0
            A[4+(f-1)*4+4,poffset+4+(f-1)*4+4] = -1.0

            A[4+(f-1)*4+3:4+(f-1)*4+4, poffset+1:poffset+2] -= Kₚ*Δ/m * d[f] * d[f]'
            A[4+(f-1)*4+3:4+(f-1)*4+4, poffset+3:poffset+4] -= Kᵥ*Δ/m * d[f] * d[f]'
            A[4+(f-1)*4+3:4+(f-1)*4+4, poffset+4+(f-1)*4+1:poffset+4+(f-1)*4+2] += Kₚ*Δ/m * d[f] * d[f]'
            A[4+(f-1)*4+3:4+(f-1)*4+4, poffset+4+(f-1)*4+3:poffset+4+(f-1)*4+4] += Kᵥ*Δ/m * d[f] * d[f]'

            l[4+(f-1)*4+3:4+(f-1)*4+4] = u[4+(f-1)*4+3:4+(f-1)*4+4] = Δ*g + Kₚ*Δ/m * d[f] * b̄[f]

            A[4+(f-1)*4+3:4+(f-1)*4+4,offset-num_faces+f] = -Δ/m * d[f]
        end
        
        S_dyn = QPN.Poly(A, l, u)
        sets[set_ind] = S_dyn
        push!(dyn_sets, set_ind)
        #Q = spzeros(total_dim, total_dim)
        #q = zeros(total_dim)
        #f = Quadratic(Q,q)
        #qp = QP(f, Dict(set_ind=>1.0), collect(offset+1:offset+4*(num_faces+1)))
        #qps[qp_ind] = qp
        #push!(level_2_progs, qp_ind)
        set_ind += 1
        #qp_ind += 1

        #for f in 1:num_faces
        Q = spzeros(total_dim, total_dim)
        Q[offset-num_faces+1:offset, offset-num_faces+1:offset] = I(num_faces)
        q = zeros(total_dim)
        f = Quadratic(Q,q)
        A_priv = spzeros(num_surfaces*num_faces, total_dim)
        l = zeros(num_surfaces*num_faces)
        for s in 1:num_surfaces
            for f in 1:num_faces
                A_priv[(s-1)*num_faces+f, offset+4+(f-1)*4+1:offset+4+(f-1)*4+2] = surface_normals[s]
                l[(s-1)*num_faces+f] = surface_nominals[s]
            end
        end

        S_priv = QPN.Poly(A_priv, l, fill(Inf, num_surfaces*num_faces))
        #vars_for_surface = collect(offset-num_faces+1:offset)
        vars_for_surface = [collect(offset-num_faces+1:offset); collect(offset+1:offset+4*(num_faces+1))]
        sets[set_ind] = S_priv
        qp = QP(f, Dict(set_ind=>1.0, set_ind-1=>1.0), vars_for_surface)
        qps[qp_ind] = qp
        push!(level_1_progs, qp_ind)
        set_ind += 1
        qp_ind += 1
    end

    #net = [level_1_progs, level_2_progs]
    net = [level_1_progs, ]
    qp_net = QPNet(qps, sets, net)
    @infiltrate

    x, Sol = solve(qp_net, x; debug=true, gen_Sol=false, high_dim=true)
   
    # setup visualization
    f = Figure()
    ax = f[1, 1] = Axis(f, aspect = DataAspect())
    xlims!(ax, -5.0, 5.0)
    ylims!(ax, 0.0, 10.0)
    p = Observable(Circle(Point(p0...), 0.1f0))
    v = Observable([Point(0.0,0.0) for i in 0:num_faces])
    verts = [x0[s:s+1] for s in 4+1:4:length(x0)]
    push!(verts, verts[1])
    v[] = verts
    lines!(ax, v, color=:red)
    poly!(ax, p, color=:blue)
    for s in 1:length(surface_normals)
        x1 = -5.0
        x2 = 5.0
        y1 = (surface_nominals[s] - surface_normals[s][1]*x1) / surface_normals[s][2]
        y2 = (surface_nominals[s] - surface_normals[s][1]*x2) / surface_normals[s][2]
        lines!(ax, [Point2f(x1,y1), Point2f(x2,y2)], color=:black)
    end
    display(f)

    pp = [x[s:s+1] for s in 1:sim_state_dim:total_dim]
    vv = [x[s+2:s+3] for s in 1:sim_state_dim:total_dim]
    λ = [x[s-num_faces+1:s] for s in sim_state_dim:sim_state_dim:total_dim]
    all_verts = [[x[t+s:t+s+1] for s in 4+1:4:length(x0)] for t in 0:sim_state_dim:T*sim_state_dim]
    for verts in all_verts
        push!(verts, verts[1])
    end

    @infiltrate
    for t = 1:length(pp)
        v[] = all_verts[t]
        p[] = Circle(Point(pp[t]...), 0.1f0)
        display(f)
        sleep(1.0)
    end

    #    pp = x[sim_state_dim+1:sim_state_dim+2]
    #    b = x[sim_state_dim+4+1:2:sim_state_dim+4+2*num_faces]
    #    verts = compute_vertices(pp, poly_normals, b)
    #    push!(verts, verts[1])
    #    v[] = verts
    #    p[] = Circle(Point(pp...), 0.1f0)
    #    display(f)
    #    x = [x[sim_state_dim+1:end]; zeros(sim_state_dim)]
    #    ProgressMeter.next!(prog, spinner="🌑🌒🌓🌔🌕🌖🌗🌘")
    #    sleep(0.001)
    #end
end


