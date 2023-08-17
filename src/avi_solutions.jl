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


"""
Assumes P has slices defined in a space described per the following ordering
[ var_vars | extra_vars | param_vars]

After permuting, slices have var_vars at locations specifed by var_inds,
param_vars are at locations specified by param_inds, and extra_vars
occupy remaining locations.
"""
function permute!(P::Poly, var_inds, param_inds)
    d = embedded_dim(P)
    dv = length(var_inds)
    dp = length(param_inds)
    extra_inds = setdiff(1:d, [var_inds ; param_inds])
    for slice in P
        a = similar(slice.a)
        a[var_inds] = slice.a[1:dv]
        a[param_inds] = slice.a[d-dp+1:d]
        a[extra_inds] = slice.a[dv+1:d-dp]
        slice.a .= a
        dropzeros!(slice.a)
    end
end

function unpermute(request::Set{Linear}, dim, var_inds, param_inds)
    isempty(request) && return request
    example = first(request)
    dv = length(var_inds)
    dp = length(param_inds)
    extra_inds = setdiff(1:dim, [var_inds; param_inds])

    appropriate_requests = Iterators.filter(request) do req
        length(req.a) == dim
    end
 
    original_request = Iterators.map(appropriate_requests) do req
        a = req.a
        a_orig = zeros(dim)
        a_orig[1:dv] = a[var_inds]
        a_orig[dv+1:dim-dp] = a[extra_inds]
        a_orig[end-dp+1:end] = a[param_inds]
        Linear(a_orig)
    end |> Set 
end

function project_and_permute(S, var_inds, param_inds)
    d = embedded_dim(S)
    dv = length(var_inds)
    dp = length(param_inds)
    projection_inds = [collect(1:dv); collect(d-dp+1:d)]
   
    piece = project(S, projection_inds)
    permute!(piece.poly, var_inds, param_inds) 
    permute!(piece.parent, var_inds, param_inds)

    return simplify(piece)
end

mutable struct LocalGAVISolutions
    gavi::GAVI
    z::Vector{Float64}
    w::Vector{Float64}
    level::Int
    subpiece_index::Int
    unexplored_Ks::Set{PolyRecipe}
    explored_Ks::Set{PolyRecipe}
    unexplored_vertices::Set{Vertex}
    explored_vertices::Set{Vertex}
    max_vertices::Int
    polys::Set{Poly}
    decision_inds::Vector{Int}
    param_inds::Vector{Int}
    permuted_request::Set{Linear}
    LocalGAVISolutions(gavi::GAVI, 
                       z::Vector{Float64}, 
                       w::Vector{Float64}, 
                       level::Int, 
                       subpiece_index::Int, 
                       decision_inds::Vector{Int}, 
                       param_inds::Vector{Int},
                       request::Set{Linear};
                       max_vertices=typemax(Int)) = begin
        n = length(z)
        m = length(w)
        permuted_request = unpermute(request, n+m, decision_inds, param_inds)
        J = comp_indices(gavi,z,w,permuted_request)
        Ks = all_Ks(J)
        @debug "There are $(length(Ks)) immediately available pieces of this solution map." 
        polys = Set{Poly}()
        explored_Ks = Set{PolyRecipe}()
        vertex_queue = Set{Vertex}()
        v = Vertex(v=[z;w])
        explored_vertices = Set((v,))
        new(gavi, z, w, level, subpiece_index, Ks, explored_Ks, vertex_queue, explored_vertices, max_vertices, polys, decision_inds, param_inds, permuted_request)
    end
end

function potential_length(ls::LocalGAVISolutions)
    length(ls.Ks) + length(ls.explored_Ks)
end

function depth(ls::LocalGAVISolutions)
    1
end
 
function get_single_solution(gavi, z, w, level, subpiece_index, decision_inds, param_inds, rng; debug=false, extra_rounds=0, permute=true, max_walk=1000.0)
    n = length(z)
    dx = length(decision_inds) + length(param_inds)
    m = length(w)

    local piece
    local x

    @infiltrate debug

    J = comp_indices(gavi,z,w)
    #K = random_K(J, rng)
    K = max_freedom_K(J)
 
    for round in 1:extra_rounds
        q = randn(rng, n)
        (; piece, reduced_inds) = local_piece(gavi,n,m,K,level,subpiece_index)
        (; A,l,u,rl,ru) = vectorize(piece)
        Aw = A[:,n+1:end]*w
        modl = OSQP.Model()
        OSQP.setup!(modl; 
                    q,
                    A=[A[:,1:n]; q'], 
                    l = [l-Aw; -max_walk], 
                    u = [u-Aw; max_walk],
                    verbose = false,
                    eps_abs=1e-8,
                    eps_rel=1e-8,
                    polish_refine_iter=10000,
                    polish=true)
        res = OSQP.solve!(modl)
        if res.info.status_val ∈ (1,2)
            #if res.info.status_polish != 1
            #    @info "Solved but not polished. Level=$level"
            #    @infiltrate level==4
            #    continue
            #end
            if !isapprox(z, res.x; atol=1e-4)
                z .= res.x
                J = comp_indices(gavi,z,w)
                K = max_freedom_K(J)
                break
            end
        end
        round == extra_rounds && @info "Did not find less restricted solution."
    end

    nv = length(decision_inds)
    reducible_inds = nv+1:n
    (; piece, reduced_inds) = local_piece(gavi,n,m,K,level,subpiece_index;reducible_inds)
    nr = setdiff(1:n, reduced_inds)
    if permute 
        permute!(piece, decision_inds, param_inds)
    end

    x = zeros(dx)
    x[decision_inds] = z[1:length(decision_inds)]
    x[param_inds] = w

    (; piece, x_opt=x, reduced_inds, z)
end


"""
Assume that J is a dict: i=>S ⊂ {1,2,3,4,5,6,7,8,9,10,11,12}
"""
function all_Ks(J)
    N = length(J)
    Ks = map(Iterators.product([J[i] for i = 1:N]...)) do assignment
        K = Dict(j=>Set(findall(x->x==j, assignment)) for j = 1:8)
    end |> Set
end

function random_K(J, rng)
    N = length(J)
    K = Dict(i=>Set{Int}() for i = 1:8)
    for (k,v) in J
        s = rand(rng, v)
        push!(K[s], k)
    end
end

function set_guide!(gavi_sols::LocalGAVISolutions, guide)
    gavi_sols.guide = guide
    for (PE, fval) in gavi_sols.polys
        gavi_sols.polys[PE] = permute_eval(guide, PE.ex, gavi_sols.decision_inds, gavi_sols.param_inds)
    end
    for (V, fval) in gavi_sols.vertex_queue
        gavi_sols.vertex_queue[V] = permute_eval(guide, V.v, gavi_sols.decision_inds, gavi_sols.param_inds)
    end
    for (K, fval) in gavi_sols.Ks
        gavi_sols.Ks[K] = permute_eval(guide, K.ex, gavi_sols.decision_inds, gavi_sols.param_inds)
    end
end
function set_guide!(::PolyUnion, guide)
end
    
function expand(gavi,z,w,K,level,subpiece_index,decision_inds,param_inds; high_dim=false)
    n = length(z)
    m = length(w)

    nv = length(decision_inds)
    reducible_inds = []
    (; piece, reduced_inds) = local_piece(gavi,n,m,K,level,subpiece_index; reducible_inds)
    if isempty(piece)
        throw(error("Piece is empty"))
    end
    
    if [z;w] ∈ piece
        slice_recipe = [z[1:nv]; fill(missing, n-nv); w]
        (; V,R,L) = get_verts(simplify(poly_slice(piece, slice_recipe)))
        vertices = [ [z[1:nv]; v; w] for v in V ]
    else
        vertices = [] 
    end
    piece = project_and_permute(piece, decision_inds, param_inds)
    (; piece, vertices)
end

function permute_eval(guide, v, decision_inds, param_inds)
    x = zeros(length(decision_inds) + length(param_inds))
    x[decision_inds] = v[1:length(decision_inds)]
    x[param_inds] = v[end-length(param_inds)+1:end]
    guide(x)
end

function Base.IteratorSize(gavi_sols::LocalGAVISolutions)
    Base.SizeUnknown()
end
function Base.eltype(gavi_sols::LocalGAVISolutions)
    Poly
end

function Base.iterate(gavi_sols::LocalGAVISolutions)
    iter_ret = iterate(gavi_sols.polys)
    if isnothing(iter_ret)
        gavi_sol_state = (; exploration_mode=true)
        return Base.iterate(gavi_sols, gavi_sol_state)
    else 
        (poly, pq_state) = iter_ret
        gavi_sol_state = (; pq_state, exploration_mode = false)
        return (poly, gavi_sol_state)
    end
end

function Base.iterate(gavi_sols::LocalGAVISolutions, state)
    if !state.exploration_mode
        ret = Base.iterate(gavi_sols.polys, state.pq_state)
        if !isnothing(ret)
            (next, pq_state) = ret
            gavi_sol_state = (; pq_state, exploration_mode = false)
            return (next, gavi_sol_state)
        end
    end
    # exploration mode (either continuing or starting)
    gavi_sol_state = (; exploration_mode = true)
    if !isempty(gavi_sols.unexplored_Ks) # if recipes available, process
        @debug "Processing recipe. Length of queue: $(length(gavi_sols.unexplored_Ks))." 
        K = pop!(gavi_sols.unexplored_Ks)
        push!(gavi_sols.explored_Ks, K)
        try
            (; piece, vertices) = expand(gavi_sols.gavi, 
                                         gavi_sols.z, 
                                         gavi_sols.w, 
                                         K, 
                                         gavi_sols.level, 
                                         gavi_sols.subpiece_index, 
                                         gavi_sols.decision_inds, 
                                         gavi_sols.param_inds)
            push!(gavi_sols.polys, piece)
            for v in vertices
                vert = Vertex(v=v)
                if vert ∉ gavi_sols.explored_vertices
                    push!(gavi_sols.unexplored_vertices, vert)
                end
            end
            return (piece, gavi_sol_state)
        catch err
            return Base.iterate(gavi_sols, gavi_sol_state)
        end
    elseif !isempty(gavi_sols.unexplored_vertices) && length(gavi_sols.explored_vertices) < gavi_sols.max_vertices # No ready-to-process Poly recipes, need to pull from available vertices
        @debug "Exploring vertex. Length of queue: $(length(gavi_sols.unexplored_vertices)). Num explored: $(length(gavi_sols.explored_vertices))"
        v = pop!(gavi_sols.unexplored_vertices)
        push!(gavi_sols.explored_vertices, v)
        J = comp_indices(gavi_sols.gavi, v.v[1:length(gavi_sols.z)], v.v[length(gavi_sols.z)+1:end], gavi_sols.permuted_request)
        Ks = all_Ks(J) |> (Ks -> setdiff(Ks, gavi_sols.explored_Ks))
        union!(gavi_sols.unexplored_Ks, Ks)
        return Base.iterate(gavi_sols, state)
    else
        return nothing
    end
end

function find_non_trivial(A,l,u,reduced_inds)
    # Non-trivial means at least one finite bound AND non-empty constraint
    non_zero_rows = rowvals(A)
    ret = [(!isinf(l[i]) || !isinf(u[i])) && i ∈ non_zero_rows for i in 1:length(l)]
end

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
function local_piece(gavi::GAVI, n, m, K, level, subpiece_index; reducible_inds=Vector{Int}(), debug=nothing)
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
            [-gavi.o[i] Inf gavi.l1[i] gavi.l1[i]]
        elseif i ∈ K[2]
            [-gavi.o[i] -gavi.o[i] gavi.l1[i] gavi.u1[i]] 
        elseif i ∈ K[3]
            [-Inf -gavi.o[i] gavi.u1[i] gavi.u1[i]]
        elseif i ∈ K[4]
            [-Inf Inf gavi.l1[i] gavi.u1[i]]
        elseif i ∈ K[5]
            [0 Inf gavi.l2[i-d1] gavi.l2[i-d1]]
        elseif i ∈ K[6]
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
    droptol!(A, 1e-8)
   
    reduced_vals = Dict{Int, Float64}()
   
    if length(reducible_inds) > 0
        while true
            further_reduced = false
            for i in 1:size(A,1)
                J = A[i,:].nzind
                reduced_inds = keys(reduced_vals)
                already_reduced = J ∩ reduced_inds
                notyet_reduced = setdiff(J, reduced_inds)
                J_reducible = notyet_reduced ∩ reducible_inds
                if isapprox(l[i], u[i]; atol=1e-6) && length(J_reducible) == 1 && notyet_reduced == J_reducible
                    j = J_reducible[1]
                    reduced_vals[j] = (u[i] - sum(A[i,k]*reduced_vals[k] for k in already_reduced; init=0.0)) / A[i,j]
                    further_reduced = true
                end
            end
            if !further_reduced
                break
            end
        end
    
        reduced_inds = keys(reduced_vals) |> collect
        reduced_vals = values(reduced_vals) |> collect
        notreduced_inds = setdiff(1:size(A,2), reduced_inds)
        remaining_reducible = Set(notreduced_inds ∩ reducible_inds)
        remaining_primary = setdiff(notreduced_inds, reducible_inds)
    
        while true
            changed = false
            for j in remaining_reducible
                con_list = A[:,j].nzind
                if !all(A[i,:].nzind ∩ remaining_reducible == A[i,:].nzind for i in con_list)
                    delete!(remaining_reducible, j)
                    changed = true
                end
            end
            if length(remaining_reducible) == 0 || !changed
                break
            end
        end
    
        r = A[:, reduced_inds]*reduced_vals
        l -= r
        u -= r
        setdiff!(notreduced_inds, remaining_reducible)
        union!(reduced_inds, remaining_reducible)
        A = A[:, notreduced_inds]
    else
        reduced_inds = []
    end

    meaningful = find_non_trivial(A,l,u,reduced_inds)
    (; piece = simplify(Poly(A[meaningful,:], l[meaningful], u[meaningful])), reduced_inds)
end

"""

TODO update this signature

Form dictionary of index sets:
r = Mz+Nw+o
J[1] = {i : lᵢ = zᵢ     , rᵢ > 0 }
J[2] = {i : lᵢ = zᵢ     , rᵢ = 0 }
J[3] = {i : lᵢ < zᵢ < uᵢ, rᵢ = 0 }
J[4] = {i :      zᵢ = uᵢ, rᵢ = 0 }
J[5] = {i :      zᵢ = uᵢ, rᵢ < 0 }
J[6] = {i : lᵢ = zᵢ = uᵢ, rᵢ free}
"""
#function comp_indices(avi::AVI, r, z, w, permuted_request=Set{Linear}(); tol=1e-4)
function comp_indices(M, N, A, B, l, u, r, z, w, permuted_request=Set{Linear}(); tol=1e-4)
    equal_bounds = isapprox.(l, u; atol=tol)
    riszero = isapprox.(r, 0; atol=tol)
    d = size(M,2) + size(N,2)
    num_requests = length(permuted_request)
    requests_granted = 0
    J = Dict{Int, Set{Int}}()
    for i = 1:length(z)

        Ji = Int[]

        a1 = -[A[i,:]; B[i,:]] # want to increase zᵢ
        a2 = -[M[i,:]; B[i,:]] # want to increase rᵢ
        try
            (l_pos, n) = lexico_positive(a1)
            a1 ./ n
        catch err
        end
        try
            (l_pos, n) = lexico_positive(a2)
            a2 ./ n
        catch err
        end
        a3 = -a2 # want to decrease rᵢ
        a4 = -a1 # want to decrease zᵢ

        for (a, j, b) in zip((a1,a2,a3,a4), (2,1,3,2), (0.0, l[i], u[i], 0.0))
            if !isinf(b) && any(a ≈ req.a for req in permuted_request)
                push!(Ji, j)
            end
        end

        if isapprox(z[i], l[i]; atol=tol) && r[i] ≥ -tol && !equal_bounds[i]
            push!(Ji,1)
        end
        if l[i]-tol ≤ z[i] ≤ u[i]+tol && riszero[i] && !equal_bounds[i]
            push!(Ji,2)
        end
        if isapprox(z[i], u[i]; atol=tol) && r[i] ≤ tol && !equal_bounds[i]
            push!(Ji,3)
        end
        if isempty(Ji)
            @assert equal_bounds[i]
            push!(Ji, 4)
        end
        J[i] = Set(Ji)
    end

    #J[1] = findall( isapprox.(z, avi.l; atol=tol) .&& r .> tol .&& .!equal_bounds)
    #J[2] = findall( isapprox.(z, avi.l; atol=tol) .&& riszero .&& .!equal_bounds)
    #J[3] = findall( (avi.l.+tol .< z .< avi.u.-tol) .&& riszero .&& .!equal_bounds)
    #J[4] = findall( isapprox.(z, avi.u; atol=tol) .&& riszero .&& .!equal_bounds)
    #J[5] = findall( isapprox.(z, avi.u; atol=tol) .&& r .< -tol .&& .!equal_bounds)
    #J[6] = findall( equal_bounds)
    #sum(length.(values(J))) < length(z) && throw(error("Z does not cleanly solve AVI"))
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
J[2] = {i : l1ᵢ = z1ᵢ      , r1ᵢ = 0 } (WEAK)
J[3] = {i : l1ᵢ < z1ᵢ < u1ᵢ, r1ᵢ = 0 }
J[4] = {i :       z1ᵢ = u1ᵢ, r1ᵢ = 0 } (WEAK)
J[5] = {i :       z1ᵢ = u1ᵢ, r1ᵢ < 0 }
J[6] = {i : l1ᵢ = z1ᵢ = u1ᵢ, r1ᵢ free}
J[7]  = {i+d1 : l2ᵢ = s2       , r2ᵢ > 0 }
J[8]  = {i+d1 : l2ᵢ = s2ᵢ      , r2ᵢ = 0 } (WEAK)
J[9]  = {i+d1 : l2ᵢ < s2ᵢ < u2ᵢ, r2ᵢ = 0 }
J[10] = {i+d1 :       s2ᵢ = u2ᵢ, r2ᵢ = 0 } (WEAK)
J[11] = {i+d1 :       s2ᵢ = u2ᵢ, r2ᵢ < 0 }
J[12] = {i+d1 : l2ᵢ = s2ᵢ = u2ᵢ, r2ᵢ free}
"""
function comp_indices(gavi::GAVI, z, w, permuted_request=Set{Linear}(); tol=1e-4)
    d1 = length(gavi.o)
    d2 = length(gavi.l2)
    @assert length(z) == d1+d2
    m = length(w)

    r1 = gavi.M*z+gavi.N*w+gavi.o
    z1 = z[1:d1]
    J1 = comp_indices(gavi.M, gavi.N, sparse(1.0I, d1,d1+d2), spzeros(d1,m), gavi.l1, gavi.u1, r1, z1, w, permuted_request; tol)
   
    M2 = [spzeros(d2, d1) sparse(I,d2,d2)]
    N2 = spzeros(d2, length(w))
    o2 = zeros(d2)
    r2 = M2*z
    s2 = gavi.A*z+gavi.B*w
    J2 = comp_indices(M2, N2, gavi.A, gavi.B, gavi.l2, gavi.u2, r2, s2, w, permuted_request; tol)

    J = Dict{Int, Set{Int}}()
    for (key, value) in J1
        J[key] = value
    end
    for (key, value) in J2
        J[key+d1] = Set(value .+ 4)
    end
    J
end
