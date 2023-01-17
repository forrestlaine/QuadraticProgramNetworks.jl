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
 
function get_single_avi_solution(avi, z, w, decision_inds, param_inds, rng; debug=false, extra_rounds=0)
    n = length(z)
    dx = length(decision_inds) + length(param_inds)
    m = length(w)

    local piece
    local x
    
    for round in 1:extra_rounds
        q = randn(rng, n)
        J = comp_indices(avi,z,w)
        K = random_K(J, rng)
        piece = local_piece(avi,n,m,K)
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
    piece = local_piece(avi,n,m,K)
    permute!(piece, decision_inds, param_inds)

    x = zeros(dx)
    x[decision_inds] = z[1:length(decision_inds)]
    x[param_inds] = w

    (; piece, x_opt=x)
end

function random_K(J, rng)
    n2 = length(J[2])
    n4 = length(J[4])
    i2 = rand(rng, Bool, n2)
    i4 = rand(rng, Bool, n4)
    K1 = Set([J[1]; [j for (i, j) in zip(i2, J[2]) if i]])
    K2 = Set([J[3]; [j for (i, j) in zip(i2, J[2]) if !i]; [j for (i, j) in zip(i4, J[4]) if !i]])
    K3 = Set([J[5]; [j for (i, j) in zip(i4, J[4]) if i]])
    K4 = Set(J[6])
    Dict(1=>K1, 2=>K2, 3=>K3, 4=>K4)
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
    piece = local_piece(avi,n,m,K)
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

function find_non_trivial(A,l,u)
    [!(isinf(l[i]) && isinf(u[i])) || all(iszero.(A[i,:])) for i in 1:length(l)]
end

"""
K[1] : Mz+Nw+o ≥ 0, z = l
K[2] : Mz+Nw+o = 0, l ≤ z ≤ u
K[3] : Mz+Nw+o ≤ 0, z = u
K[4] : Mz+Nw+o free, l = z = u
"""
function local_piece(avi, n, m, K)
    A = [avi.M avi.N;
         I(n) spzeros(n,m)]

    bounds = mapreduce(vcat, 1:n) do i
        if i ∈ K[1]
            [-avi.o[i] Inf avi.l[i] avi.l[i]]
        elseif i ∈ K[2]
            [-avi.o[i] -avi.o[i] avi.l[i] avi.u[i]] 
        elseif i ∈ K[3]
            [-Inf -avi.o[i] avi.u[i] avi.u[i]]
        else
            [-Inf Inf avi.l[i] avi.u[i]]
        end
    end
    l = [bounds[:,1]; bounds[:,3]]
    u = [bounds[:,2]; bounds[:,4]]
    inds = l.>u
    l[inds] = u[inds]
    meaningful = find_non_trivial(A,l,u)
    Poly(A[meaningful,:], l[meaningful], u[meaningful])
end

"""
Form dictionary of index sets:
J[1] = {i : lᵢ = zᵢ     , (Mz+Nw+o)ᵢ > 0 }
J[2] = {i : lᵢ = zᵢ     , (Mz+Nw+o)ᵢ = 0 }
J[3] = {i : lᵢ < zᵢ < uᵢ, (Mz+Nw+o)ᵢ = 0 }
J[4] = {i :      zᵢ = uᵢ, (Mz+Nw+o)ᵢ = 0 }
J[5] = {i :      zᵢ = uᵢ, (Mz+Nw+o)ᵢ < 0 }
J[6] = {i : lᵢ = zᵢ = uᵢ, (Mz+Nw+o)ᵢ = 0 }
"""
function comp_indices(avi, z, w; tol=1e-4)
    J = Dict{Int, Vector{Int}}()
    r = avi.M*z+avi.N*w+avi.o
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
