PolyRecipe = Dict{Int, Set{Int}}

"""
A polyhedral region along with an accompanying 
exemplar point ex.
ex may be a vertex, centroid, or other, and should
not be assumed to have any specific properties other than
that ex ∈ poly.
"""
struct PolyExemplar
    poly::Poly
    ex::Vector{Float64}
end

struct RecipeExemplar
    recipe::PolyRecipe
    ex::Vector{Float64}
end

"""
Simple wrapper around Vector{Float64} to enable numerical 'isequal'.
"""
Base.@kwdef struct Vertex
    v::Vector{Float64}
    digits::Int=5
end
function Base.isequal(v1::Vertex, v2::Vertex)
    isequal(round.(v1.v; digits=v1.digits), round.(v2.v;digits=v2.digits))
end
function Base.hash(v::Vertex, h::UInt)
    hash(round.(v.v; digits=v.digits))
end

function permute!(P, var_inds, param_inds)
    d = embedded_dim(P)
    dv = length(var_inds)
    dp = length(param_inds)
    for slice in P
        a = similar(slice.a)
        a[var_inds] = slice.a[1:dv]
        a[param_inds] = slice.a[d-dp+1:d]
        a[dv+dp+1:end] = slice.a[dv+1:d-dp]
        slice.a .= a
        dropzeros!(slice.a)
    end
end

function project_and_permute(S, var_inds, param_inds)
    d = embedded_dim(S)
    dv = length(var_inds)
    dp = length(param_inds)
    projection_inds = [collect(1:dv); collect(d-dp+1:d)]
    
    piece = project(S, projection_inds)
    permute!(piece, var_inds, param_inds) 

    return simplify(piece)
end

mutable struct LocalAVISolutions
    avi::AVI
    z::Vector{Float64}
    w::Vector{Float64}
    guide::Function
    vertex_queue::PriorityQueue{Vertex, Float64}
    Ks::PriorityQueue{RecipeExemplar, Float64}
    explored_vertices::Set{Vertex}
    explored_Ks::Set{PolyRecipe}
    polys::PriorityQueue{PolyExemplar, Float64}
    decision_inds::Vector{Int}
    param_inds::Vector{Int}
    LocalAVISolutions(avi::AVI, z::Vector{Float64}, w::Vector{Float64}, decision_inds::Vector{Int}, param_inds::Vector{Int}) = begin
        n = length(z)
        m = length(w)
        J = comp_indices(avi,z,w)
        J2 = Set(J[2])
        J4 = Set(J[4])
        Ks = map(enumerate(Iterators.product(powerset(J[2]), powerset(J[4])))) do (e,(S2, S4))
            C2 = setdiff(J2, Set(S2)) |> collect
            C4 = setdiff(J4, Set(S4)) |> collect
            Dict(1=>Set([J[1];C2]), 2=>Set([J[3];S2;S4]),3=>Set([J[5];C4]),4=>Set(J[6]))
        end |> Set
        K = pop!(Ks)
        (; piece, exemplar, vertices) = expand(avi,z,w,K,decision_inds,param_inds)
        polys = PriorityQueue{PolyExemplar, Float64}()
        vertex_queue = PriorityQueue{Vertex, Float64}()
        enqueue!(polys, PolyExemplar(piece, exemplar), Inf)
        for v in vertices
            enqueue!(vertex_queue, Vertex(v=v), Inf)
        end
        queued_Ks = PriorityQueue{RecipeExemplar, Float64}()
        for KK in Ks
            queued_Ks[RecipeExemplar(KK, [z;w])] = Inf
        end
        explored_Ks = Set{PolyRecipe}((K,))
        explored_vertices = Set{Vertex}()
        guide = (x->Inf)
        new(avi, z, w, guide, vertex_queue, queued_Ks, explored_vertices, explored_Ks, polys, decision_inds, param_inds)
    end
end
 
function get_single_avi_solution(avi, z, w, decision_inds, param_inds, rng; debug=false, extra_rounds=0, permute=true)
    n = length(z)
    dx = length(decision_inds) + length(param_inds)
    m = length(w)

    local piece
    local x
    
    for round in 1:extra_rounds
        q = randn(rng, n)
        J = comp_indices(avi,z,w)
        K = random_K(J, rng)
        (; piece, reduced_inds) = local_piece(avi,n,m,K)
        (A,l,u,rl,ru) = vectorize(piece)
        Aw = A[:,n+1:end]*w
        mod = OSQP.Model()
        OSQP.setup!(mod; 
                    q,
                    A=[A[:,1:n]; q'], 
                    l = [l-Aw; -10.0], 
                    u = [u-Aw; 10.0],
                    verbose = false,
                    eps_abs=1e-8,
                    eps_rel=1e-8,
                    polish=true)
        res = OSQP.solve!(mod)
        if res.info.status_val == 1
            z = res.x
        end
    end

    J = comp_indices(avi,z,w)
    K = random_K(J, rng)

    nv = length(decision_inds)
    np = length(param_inds)
    nd = n - nv - np
    reducible_inds = nv+1:n
    (; piece, reduced_inds) = local_piece(avi,n,m,K; reducible_inds)
    if permute 
        permute!(piece, decision_inds, param_inds)
    end

    x = zeros(dx)
    x[decision_inds] = z[1:length(decision_inds)]
    x[param_inds] = w

    (; piece, x_opt=x, reduced_inds)
end


#"""
#    [B1 B2][x; y] = b
#l ≤ [C1 C2][x; y] ≤ u
#
#size(B2) = m × n
#if m ≥ n && rank(B2) == m
#    y can be eliminated trivially
#if m ≥ n && m > rank(B2) ≥ n
#    y can be partially eliminated trivially (I THINK)
#if m < n
#    need to think this through...
#
#"""
#function partial_project(piece, dx; tol=1e-4, sval_tol=1e-6, sp_tol=1e-8)
#    (A,l,u,_,_) = vectorize(piece)
#    eq_inds = collect(1:length(l))[isapprox.(l, u; atol=tol)]
#    iq_inds = setdiff(1:length(l), eq_inds) |> sort
#    Ae = A[eq_inds,:]
#    Ai = A[iq_inds,:]
#    b = u[eq_inds]
#    A1 = Ae[:,1:dx]
#    A2 = Ae[:,dx+1:end]
#    r = rank(A2)
#    if r == size(A2,2) # linearly independent columns
#        Atmp = A2'A2
#        pinv = mapreduce(hcat, eachrow(A2)) do a
#            sparse(Atmp \ collect(a))
#        end
#        resid = A2*pinv - I
#        droptol!(resid, sp_tol)
#        A_eq_new = resid * A1
#        b_eq_new = resid * b
#
#        A_ineq_new = Ai[:,1:dx] - Ai[:,dx+1:end] * pinv*A1
#        tmp = Ai[:,dx+1:end]*pinv*b
#        u_ineq_new = u[iq_inds]-tmp
#        l_ineq_new = l[iq_inds]-tmp
#        p = Poly([A_eq_new; A_ineq_new], [b_eq_new; l_ineq_new], [b_eq_new; u_ineq_new])
#        @infiltrate
#        return p
#    else
#        @infiltrate
#        return piece
#    end
#end

"""
K[1] : J1 ∪ J2a
K[2] : J2b ∪ J3 ∪ J4a
K[3] : J4b ∪ J5
K[4] : J6

IF J[7-12] exist:

K[5] : J7 ∪ J8a
K[6] : J8b ∪ J9 ∪ J10a
K[7] : J10b ∪ J11
K[8] : J12

Jia/Jib is random partitioning of Ji
"""
function random_K(J, rng)
    n2 = length(J[2])
    n4 = length(J[4])
    i2 = rand(rng, Bool, n2)
    i4 = rand(rng, Bool, n4)
    K1 = Set([J[1]; [j for (i, j) in zip(i2, J[2]) if i]])
    K2 = Set([J[3]; [j for (i, j) in zip(i2, J[2]) if !i]; [j for (i, j) in zip(i4, J[4]) if !i]])
    K3 = Set([J[5]; [j for (i, j) in zip(i4, J[4]) if i]])
    K4 = Set(J[6])
    K = Dict(1=>K1, 2=>K2, 3=>K3, 4=>K4)

    if 7 ∈ keys(J) # J are indices corresponding to GAVI solution
        n8 = length(J[8])
        n10 = length(J[10])
        i8 = rand(rng, Bool, n8)
        i10 = rand(rng, Bool, n10)
        K[5] = Set([J[7]; [j for (i, j) in zip(i8, J[8]) if i]])
        K[6] = Set([J[9]; [j for (i, j) in zip(i8, J[8]) if !i]; [j for (i, j) in zip(i10, J[10]) if !i]])
        K[7] = Set([J[11]; [j for (i, j) in zip(i10, J[10]) if i]])
        K[8] = Set(J[12])
    end
    K
end

function set_guide!(avi_sols::LocalAVISolutions, guide)
    avi_sols.guide = guide
    for (PE, fval) in avi_sols.polys
        avi_sols.polys[PE] = permute_eval(guide, PE.ex, avi_sols.decision_inds, avi_sols.param_inds)
    end
    for (V, fval) in avi_sols.vertex_queue
        avi_sols.vertex_queue[V] = permute_eval(guide, V.v, avi_sols.decision_inds, avi_sols.param_inds)
    end
    for (K, fval) in avi_sols.Ks
        avi_sols.Ks[K] = permute_eval(guide, K.ex, avi_sols.decision_inds, avi_sols.param_inds)
    end
end
function set_guide!(::Vector{Poly}, guide)
end

function expand(avi,z,w,K,decision_inds,param_inds; high_dim=false)
    n = length(z)
    m = length(w)
    piece = local_piece(avi,n,m,K,decision_inds,param_inds)
    (; V,R,L) = get_verts(simplify(poly_slice(piece, [fill(missing, n); w])))
    vertices = [ [v;w] for v in V]
    rays = [ [r.a;zero(w)] for r in R]
    lines = [ [l.a;zero(w)] for l in L]
    avg_vertex = sum(vertices) / length(vertices)
    exemplar = avg_vertex + sum(rays; init=zeros(n+m)) + sum(lines; init=zeros(n+m))
    (; piece=project_and_permute(piece, decision_inds, param_inds), exemplar, vertices)
end

function permute_eval(guide, v, decision_inds, param_inds)
    x = zeros(length(decision_inds) + length(param_inds))
    x[decision_inds] = v[1:length(decision_inds)]
    x[param_inds] = v[end-length(param_inds)+1:end]
    guide(x)
end

function Base.IteratorSize(avi_sols::LocalAVISolutions)
    Base.SizeUnknown()
end
function Base.eltype(avi_sols::LocalAVISolutions)
    Poly
end

function Base.iterate(avi_sols::LocalAVISolutions)
    (next, pq_state) = Base.iterate(avi_sols.polys)
    @assert !isnothing(next) # avi_sols should always have at least one piece, including at initialization
    avi_sol_state = (; pq_state, exploration_mode = false)
    if !isempty(next.first.poly)
        return (next.first.poly, avi_sol_state)
    else
        return Base.iterate(avi_sols, avi_sol_state)
    end
end

function Base.iterate(avi_sols::LocalAVISolutions, state)
    if !state.exploration_mode
        ret = Base.iterate(avi_sols.polys, state.pq_state)
        if !isnothing(ret)
            (next, pq_state) = ret
            avi_sol_state = (; pq_state, exploration_mode = false)
            if !isempty(next.first.poly)
                return (next.first.poly, avi_sol_state)
            else
                return Base.iterate(avi_sols, avi_sol_state)
            end
        end
    end
    # exploration mode (either continuing or starting)
    if !isempty(avi_sols.Ks) # if recipes available, process
        K = dequeue!(avi_sols.Ks)
        (; piece, exemplar, vertices) = expand(avi_sols.avi, avi_sols.z, avi_sols.w, K.recipe, avi_sols.decision_inds, avi_sols.param_inds)
        fval = permute_eval(avi_sols.guide, exemplar, avi_sols.decision_inds, avi_sols.param_inds)
        enqueue!(avi_sols.polys, PolyExemplar(piece, exemplar), fval) 
        push!(avi_sols.explored_Ks, K.recipe)
        for v in vertices
            vert = Vertex(v=v)
            if vert ∉ avi_sols.explored_vertices
                fval = permute_eval(avi_sols.guide, v, avi_sols.decision_inds, avi_sols.param_inds)
                avi_sols.vertex_queue[vert] = fval
            end
        end
        avi_sol_state = (; exploration_mode = true)
        if !isempty(piece)
            return (piece, avi_sol_state)
        else
            return Base.iterate(avi_sols, avi_sol_state)
        end
    elseif !isempty(avi_sols.vertex_queue) # No ready-to-process Poly recipes, need to pull from available vertices
        v = dequeue!(avi_sols.vertex_queue)
        push!(avi_sols.explored_vertices, v)
        J = comp_indices(avi_sols.avi, v.v[1:length(avi_sols.z)], v.v[length(avi_sols.z)+1:end]) 
        J2 = Set(J[2])
        J4 = Set(J[4])
        Ks = map(enumerate(Iterators.product(powerset(J[2]), powerset(J[4])))) do (e,(S2, S4))
            C2 = setdiff(J2, Set(S2)) |> collect
            C4 = setdiff(J4, Set(S4)) |> collect
            Dict(1=>Set([J[1];C2]), 2=>Set([J[3];S2;S4]),3=>Set([J[5];C4]),4=>Set(J[6]))
        end |> Set |> (Ks -> setdiff(Ks,avi_sols.explored_Ks))
        for K in Ks
            fval = permute_eval(avi_sols.guide, v.v, avi_sols.decision_inds, avi_sols.param_inds)
            avi_sols.Ks[RecipeExemplar(K, v.v)] = fval
        end
        #union!(avi_sols.Ks, Ks)
        return Base.iterate(avi_sols, state)
    else
        return nothing
    end
end

function find_non_trivial(A,l,u,reduced_inds)
    # Non-trivial means at least one finite bound AND non-empty constraint
    non_zero_rows = rowvals(A)
    ret = [(!isinf(l[i]) || !isinf(u[i])) && i ∈ non_zero_rows for i in 1:length(l)]
end

#"""
#K[1] : Mz+Nw+o ≥ 0, z = l
#K[2] : Mz+Nw+o = 0, l ≤ z ≤ u
#K[3] : Mz+Nw+o ≤ 0, z = u
#K[4] : Mz+Nw+o free, l = z = u
#"""
#function local_piece(avi::AVI, n, m, K; reducible_inds=Vector{Int}())
#    A = [avi.M avi.N;
#         I(n) spzeros(n,m)]
#
#    lo_reduced = []
#    up_reduced = []
#    bounds = mapreduce(vcat, 1:n) do i
#        if i ∈ K[1]
#            i ∈ reducible_inds && push!(lo_reduced, i)
#            [-avi.o[i] Inf avi.l[i] avi.l[i]]
#        elseif i ∈ K[2]
#            [-avi.o[i] -avi.o[i] avi.l[i] avi.u[i]] 
#        elseif i ∈ K[3]
#            i ∈ reducible_inds && push!(up_reduced, i)
#            [-Inf -avi.o[i] avi.u[i] avi.u[i]]
#        else
#            i ∈ reducible_inds && push!(lo_reduced, i)
#            [-Inf Inf avi.l[i] avi.u[i]]
#        end
#    end
#    l = [bounds[:,1]; bounds[:,3]]
#    u = [bounds[:,2]; bounds[:,4]]
#    noisy_inds = l.>u
#    l[noisy_inds] = u[noisy_inds]
#
#    @infiltrate
#
#    reduced_inds = [lo_reduced; up_reduced]
#    notreduced_inds = setdiff(1:size(A,2), reduced_inds)
#    Al = A[:,lo_reduced]
#    Au = A[:,up_reduced]
#    A = A[:,notreduced_inds]
#
#    reduced_contributions = Al * bounds[lo_reduced,3] + Au * bounds[up_reduced,4]
#    l -= reduced_contributions
#    u -= reduced_contributions
#    
#    meaningful = find_non_trivial(A,l,u, reduced_inds)
#    (; piece = Poly(A[meaningful,:], l[meaningful], u[meaningful]), reduced_inds)
#end
"""
K[1] : Mz+Nw+o ≥ 0, z = l
K[2] : Mz+Nw+o = 0, l ≤ z ≤ u
K[3] : Mz+Nw+o ≤ 0, z = u
K[4] : Mz+Nw+o free, l = z = u
K[5] : z ≥ 0, Az+Bw = l
K[6] : z = 0, l ≤ Az+Bw ≤ u
K[7] : z ≤ 0, Az+Bw = u
K[8] : z free, l = Az+Bw = u
"""
function local_piece(gavi::GAVI, n, m, K; reducible_inds=Vector{Int}())
    d1 = length(gavi.l1)
    d2 = length(gavi.l2)
    I1 = [sparse(I,d1,d1) spzeros(d1,d2)]
    I2 = [spzeros(d2,d1) sparse(I,d2,d2)]
    A = [gavi.M gavi.N;
         I2 spzeros(d2,m);
         I1 spzeros(d1,m);
         gavi.A gavi.B]

    lo_reduced = []
    up_reduced = []
    zero_reduced = []
    bounds = mapreduce(vcat, 1:n) do i
        if i ∈ K[1]
            i ∈ reducible_inds && push!(lo_reduced, i)
            [-gavi.o[i] Inf gavi.l1[i] gavi.l1[i]]
        elseif i ∈ K[2]
            [-gavi.o[i] -gavi.o[i] gavi.l1[i] gavi.u1[i]] 
        elseif i ∈ K[3]
            i ∈ reducible_inds && push!(up_reduced, i)
            [-Inf -gavi.o[i] gavi.u1[i] gavi.u1[i]]
        elseif i ∈ K[4]
            i ∈ reducible_inds && push!(lo_reduced, i)
            [-Inf Inf gavi.l1[i] gavi.u1[i]]
        elseif i ∈ K[5]
            #TODO CHECK INDEX OF THESE
            [0 Inf gavi.l2[i-d1] gavi.l2[i-d1]]
        elseif i ∈ K[6]
            i ∈ reducible_inds && push!(zero_reduced, i)
            [0.0 0 gavi.l2[i-d1] gavi.u2[i-d1]]
        elseif i ∈ K[7]
            [-Inf 0 gavi.u2[i-d1] gavi.u2[i-d1]]
        elseif i ∈ K[8]
            [-Inf Inf gavi.l2[i-d1] gavi.u2[i-d1]]
        else
            @infiltrate
        end
    end
    l = [bounds[:,1]; bounds[:,3]]
    u = [bounds[:,2]; bounds[:,4]]
    noisy_inds = l.>u
    l[noisy_inds] = u[noisy_inds]

    reduced_inds = [lo_reduced; up_reduced; zero_reduced]
    notreduced_inds = setdiff(1:size(A,2), reduced_inds)
    Al = A[:,lo_reduced]
    Au = A[:,up_reduced]
    A = A[:,notreduced_inds]

    reduced_contributions = Al * bounds[lo_reduced,3] + Au * bounds[up_reduced,4] # zero_reduced inds don't contribute
    l -= reduced_contributions
    u -= reduced_contributions
   
    droptol!(A, 1e-8)
    meaningful = find_non_trivial(A,l,u, reduced_inds)
    (; piece = Poly(A[meaningful,:], l[meaningful], u[meaningful]), reduced_inds)
end

"""
Form dictionary of index sets:
r = Mz+Nw+o
J[1] = {i : lᵢ = zᵢ     , rᵢ > 0 }
J[2] = {i : lᵢ = zᵢ     , rᵢ = 0 }
J[3] = {i : lᵢ < zᵢ < uᵢ, rᵢ = 0 }
J[4] = {i :      zᵢ = uᵢ, rᵢ = 0 }
J[5] = {i :      zᵢ = uᵢ, rᵢ < 0 }
J[6] = {i : lᵢ = zᵢ = uᵢ, rᵢ = 0 }
"""
function comp_indices(avi::AVI, r, z, w; tol=1e-4)
    J = Dict{Int, Vector{Int}}()
    equal_bounds = isapprox.(avi.l, avi.u; atol=tol)
    riszero = isapprox.(r, 0; atol=tol)
    J[1] = findall( isapprox.(z, avi.l; atol=tol) .&& r .> tol )
    J[2] = findall( isapprox.(z, avi.l; atol=tol) .&& riszero .&& .!equal_bounds)
    J[3] = findall( (avi.l.+tol .< z .< avi.u.-tol) .&& riszero )
    J[4] = findall( isapprox.(z, avi.u; atol=tol) .&& riszero .&& .!equal_bounds)
    J[5] = findall( isapprox.(z, avi.u; atol=tol) .&& r .< -tol )
    J[6] = findall( equal_bounds .&& riszero )
    @infiltrate sum(length.(values(J))) != length(z)
    return J
end
function comp_indices(avi::AVI, z, w; tol=1e-4)
    r = avi.M*z+avi.N*w+avi.o
    comp_indices(avi, r, z, w; tol)
end

"""
Form dictionary of index sets:
r1 = Mz+Nw+o
d1 = length(o)
r2 = z2
s2 = Az+Bw
J[1] = {i : l1ᵢ = z1ᵢ      , r1ᵢ > 0 }
J[2] = {i : l1ᵢ = z1ᵢ      , r1ᵢ = 0 }
J[3] = {i : l1ᵢ < z1ᵢ < u1ᵢ, r1ᵢ = 0 }
J[4] = {i :       z1ᵢ = u1ᵢ, r1ᵢ = 0 }
J[5] = {i :       z1ᵢ = u1ᵢ, r1ᵢ < 0 }
J[6] = {i : l1ᵢ = z1ᵢ = u1ᵢ, r1ᵢ = 0 }
J[7]  = {i+d1 : l2ᵢ = s2       , r2ᵢ > 0}
J[8]  = {i+d1 : l2ᵢ = s2ᵢ      , r2ᵢ = 0}
J[9]  = {i+d1 : l2ᵢ < s2ᵢ < u2ᵢ, r2ᵢ = 0}
J[10] = {i+d1 :       s2ᵢ = u2ᵢ, r2ᵢ = 0}
J[11] = {i+d1 :       s2ᵢ = u2ᵢ, r2ᵢ < 0}
J[12] = {i+d1 : l2ᵢ = s2ᵢ = u2ᵢ, r2ᵢ = 0}
"""
function comp_indices(gavi::GAVI, z, w; tol=1e-4)
    avi1 = AVI(gavi.M, gavi.N, gavi.o, gavi.l1, gavi.u1)
    r1 = gavi.M*z+gavi.N*w+gavi.o
    z1 = z[1:length(gavi.o)]
    J1 = comp_indices(avi1, r1, z1, w; tol)
   
    d1 = length(gavi.o)
    d2 = length(gavi.l2)
    M2 = [spzeros(d2, d1) sparse(I,d2,d2)]
    N2 = spzeros(d2, length(w))
    o2 = zeros(d2)
    avi2 = AVI(M2, N2, o2, gavi.l2, gavi.u2)
    r2 = M2*z
    s2 = gavi.A*z+gavi.B*w
    J2 = comp_indices(avi2, r2, s2, w; tol)

    J = Dict{Int, Vector{Int}}()
    for (key,value) in J1
        J[key] = value
    end
    for (key, value) in J2
        J[key+6] = value .+ d1
    end
    for i = 1:length(z)
        if !any(i ∈ J[j] for j in 1:12)
            @infiltrate
        end
    end
    J
end
