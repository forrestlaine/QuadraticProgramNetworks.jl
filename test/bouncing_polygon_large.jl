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
    T = 10   # number of simulation steps
    Δ = 0.10  # simulation timestep
    p0 = [0.0; 4] # initial configuration
    v0 = [0.0; 2.0]  # initial velocity
    g = [0.0, -0.1]  # gravity vector
    Kₚ = 50.0 
    Kᵥ = 1.0
    M = 5.0 # centroid mass
    m = 0.5 # face mass

    num_faces = 5
    poly_normals = [[cos(θ), sin(θ)] for θ in collect(1:num_faces)*(2*pi)/num_faces]
    poly_nominals = ones(Float64, num_faces)
    surface_normals = [[1.0, 3.0], [-1.0, 1.0], [1.0, -1.0], [-1.0, -3.0]]
    surface_nominals = [3.5, 0, -8.5, -20.0]


    # poly_normal'*(r) ≤ poly_nominal
    # surface_normal'*(p+r) ≥ surface_nominal

    num_faces = length(poly_nominals)
    num_surfaces = length(surface_nominals)
    
    # setup looks like this:
    #  □    
    # \  /  
    #  \/   
  
    # Centroid: 4
    # Faces:
    #   b ∈ ℝ, d ∈ ℝ, λ ∈ ℝ
    # Surfaces:
    #   r ∈ ℝ²

    # vars: [ p0 v0 bd0 | λ1 r1 p1 v1 bd1 | λ2 r2 p2 v2 bd2 | ... ]
    sim_state_dim = num_surfaces * 2 + 4 + num_faces * 3
    total_dim = 4 + num_faces*2 + T*sim_state_dim # add initial values into state for convenience

    qps = Dict{Int, QP}()
    sets =  Dict{Int, QPN.Poly}()
    level_1_progs = Set{Int}()
    level_2_progs = Set{Int}()
    
    qp_ind = 1
    set_ind = 1

    initial_offset = 4+num_faces*2

    for t ∈ 1:T
        # Dynamics equations
        # equations for pt, vt, bdt
        a = zeros(4+num_faces*2, total_dim)
        l = zeros(4+num_faces*2)
        u = zeros(4+num_faces*2)

        offset = initial_offset+sim_state_dim * (t-1)
        poffset = offset - sim_state_dim
        surf_off = 2*num_surfaces+num_faces

        a[1, offset+surf_off+1] = 1.0
        a[1, offset+surf_off+3] = -Δ
        a[1, poffset+surf_off+1] = -1.0 
        a[2, offset+surf_off+2] = 1.0
        a[2, offset+surf_off+4] = -Δ
        a[2, poffset+surf_off+2] = -1.0
        a[3, offset+surf_off+3] = 1.0
        a[3, poffset+surf_off+3] = -1.0
        a[4, offset+surf_off+4] = 1.0 
        a[4, poffset+surf_off+4] = -1.0

        l[3] = u[3] = Δ*g[1]
        l[4] = u[4] = Δ*g[2]
        for f ∈ 1:num_faces
            a[3,poffset+surf_off+4+(f-1)*2+1] = -Δ*Kₚ*poly_normals[f][1]/M
            l[3] -= Δ*Kₚ*poly_nominals[f]*poly_normals[f][1]/M
            a[3,poffset+surf_off+4+(f-1)*2+2] = -Δ*Kᵥ*poly_normals[f][1]/M

            a[4,poffset+surf_off+4+(f-1)*2+1] = -Δ*Kₚ*poly_normals[f][2]/M
            l[4] -= Δ*Kₚ*poly_nominals[f]*poly_normals[f][2]/M
            a[4,poffset+surf_off+4+(f-1)*2+2] = -Δ*Kᵥ*poly_normals[f][2]/M
 
            a[4+(f-1)*2+1,offset+surf_off+4+(f-1)*2+1] = 1.0
            a[4+(f-1)*2+1,poffset+surf_off+4+(f-1)*2+1] = -1.0
            a[4+(f-1)*2+1,offset+surf_off+4+(f-1)*2+2] = -Δ
            a[4+f*2, offset+surf_off+4+(f-1)*2+2] = 1.0
            a[4+f*2, poffset+surf_off+4+(f-1)*2+2] = -1.0
            a[4+f*2, poffset+surf_off+4+(f-1)*2+1] += +Δ*Kₚ/m
            a[4+f*2, poffset+surf_off+4+(f-1)*2+2] += +Δ*Kᵥ/m
            l[4+f*2] = u[4+f*2] = Δ*(Kₚ*poly_nominals[f]/m + g'*poly_normals[f])
            a[4+f*2, offset+f] = Δ/m
            #for s ∈ 1:num_surfaces 
            #    a[4+f*2, offset+(s-1)*4+1] = -Δ*poly_normals[f][1]/m
            #    a[4+f*2, offset+(s-1)*4+2] = -Δ*poly_normals[f][2]/m
            #end
            #d - d̄ - Δ * (-Kp (b-nom) - Kv (d) + λ'*n) / m  = Δ
        end
        S_dyn = QPN.Poly(a, l, u)
        sets[set_ind] = S_dyn
        Q = spzeros(total_dim, total_dim)
        q = zeros(total_dim)
        f = Quadratic(Q,q)
        qp = QP(f, Dict(set_ind=>1.0), collect(offset+surf_off+1:offset+surf_off+4+num_faces*2)) # p, v, bd. No objective just satisfy dynamics
        qps[qp_ind] = qp
        push!(level_2_progs, qp_ind)
        set_ind += 1
        qp_ind += 1

        #for f in 1:num_faces
        Q = spzeros(total_dim, total_dim)
        Q[offset+1:offset+num_faces, offset+1:offset+num_faces] = I(num_faces)
        q = zeros(total_dim)
        f = Quadratic(Q,q)
        a_priv = zeros(num_surfaces, total_dim)
        for s in 1:num_surfaces
            a_priv[s,offset+num_faces+(s-1)*2+1] = a_priv[s,offset+surf_off+1] = surface_normals[s][1]
            a_priv[s,offset+num_faces+(s-1)*2+2] = a_priv[s,offset+surf_off+2] = surface_normals[s][2]
        end

        S_priv = QPN.Poly(a_priv, surface_nominals, fill(Inf, num_surfaces))
        vars_for_surface = collect(offset+1:offset+num_faces)
        sets[set_ind] = S_priv
        qp = QP(f, Dict(set_ind=>1.0), vars_for_surface)
        qps[qp_ind] = qp
        push!(level_1_progs, qp_ind)
        set_ind += 1
        qp_ind += 1
        
        for s in 1:num_surfaces
            Q = spzeros(total_dim, total_dim)
            q = zeros(total_dim)
            q[offset+num_faces+(s-1)*2+1] = surface_normals[s][1]
            q[offset+num_faces+(s-1)*2+2] = surface_normals[s][2]
            f = Quadratic(Q,q)
            a_priv = zeros(num_faces, total_dim)
            u = zeros(num_faces)
            l = fill(-Inf, num_faces)
            for f in 1:num_faces
                a_priv[f,offset+num_faces+(s-1)*2+1] = poly_normals[f][1]
                a_priv[f,offset+num_faces+(s-1)*2+2] = poly_normals[f][2]
                a_priv[f,offset+surf_off+4+(f-1)*2+1] = -1.0
            end
            S_priv = QPN.Poly(a_priv, l, u)
            sets[set_ind] = S_priv
            vars_for_r = [offset+num_faces+(s-1)*2+1, offset+num_faces+(s-1)*2+2]
            qp = QP(f, Dict(set_ind=>1.0), vars_for_r)
            qps[qp_ind] = qp
            push!(level_2_progs, qp_ind)
            set_ind += 1
            qp_ind += 1
        end
    end

    net = [level_1_progs, level_2_progs]
    qp_net = QPNet(qps, sets, net)
    
    x = [p0; v0; reduce(vcat, ([nom; 0] for nom in poly_nominals)); zeros(T*sim_state_dim)]
    
    #ps = [p0,]
    #bs = [poly_nominals,]
    #
    #f = Figure()
    #ax = f[1, 1] = Axis(f, aspect = DataAspect())
    #xlims!(ax, -5.0, 5.0)
    #ylims!(ax, 0.0, 10.0)
    #p = Observable(Circle(Point(p0...), 0.1f0))
    #v = Observable([Point(0.0,0.0) for i in 0:length(poly_normals)])
    #verts = compute_vertices(p0, poly_normals, poly_nominals)
    #push!(verts, verts[1])
    #v[] = verts
    #lines!(ax, v, color=:red)
    #poly!(ax, p, color=:blue)
    #for s in 1:length(surface_normals)
    #    x1 = -5.0
    #    x2 = 5.0
    #    y1 = (surface_nominals[s] - surface_normals[s][1]*x1) / surface_normals[s][2]
    #    y2 = (surface_nominals[s] - surface_normals[s][1]*x2) / surface_normals[s][2]
    #    lines!(ax, [Point2f(x1,y1), Point2f(x2,y2)], color=:black)
    #end
    #display(f)

    ##for t in 1:TT
    #prog = ProgressUnknown("Simulating:", spinner=true)
    #while true
    x, Sol = solve(qp_net, x; debug=true, gen_Sol=false)
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


