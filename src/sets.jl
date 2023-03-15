Relation = Union{typeof(<), typeof(≤)}

function Random.rand(rng::AbstractRNG, ::Random.SamplerType{Relation})
    r = rand(rng, Bool)
    r ? ( < ) : ( ≤ )
end

function complement(r::Relation)
    (r == ≤) ? ( < ) : ( ≤ )
end
function Base.adjoint(r::Relation)
    r
end

"""
True iff x is lexico positive.
"""
function lexico_positive(x::Vector{Float64})
    return first(x) ≥ 0 
end
function lexico_positive(x::SparseVector{Float64, Int64})
    (I, V) = findnz(x)
    t = argmin(I)
    return V[t] ≥ 0
end

"""
A normalized slice.

S := { x : l ⋈ₗ a'x ⋈ᵤ u }

a ∈ ℝⁿ is required to be of unit norm and lexico-positive.

"""
struct Slice
    a::SparseVector{Float64, Int64}
    l::Float64
    u::Float64
    rl::Relation
    ru::Relation
    Slice(a,l,u,rl,ru; tol=1e-6) = begin
        n = norm(a)
        if (n ≤ tol)
            new(spzeros(length(a)), l, u, rl, ru)
        else
            n = 1.0
            if lexico_positive(a) || true
                new(sparse(a/n), l/n, u/n, rl, ru)
            else
                new(sparse(-a/n), -u/n, -l/n, ru, rl)
            end
        end
    end
    Slice(a,l,u) = Slice(a,l,u,≤,≤)
end

"""
return the closure of the slice S
"""
function closure(s::Slice)
    Slice(copy(s.a), s.l, s.u, ≤, ≤)     
end

"""
To be used for hashing purposes -- proper defn of Polyhedral regions
"""
function Base.isequal(s1::Slice, s2::Slice) 
    (round.(s1.a; digits=5)==round.(s2.a; digits=5)) &&
        (round(s1.l; digits=5)==round(s2.l; digits=5)) &&
        (round(s1.u; digits=5) == round(s2.u; digits=5)) &&
        s1.rl==s2.rl && s1.ru == s2.ru
end
function Base.hash(s::Slice, h::UInt)
    hash((round.(s.a; digits=5), round(s.l; digits=5), round(s.u; digits=5), s.rl, s.ru), h)
end

"""
A not-necessarily closed polyhedron. 
Implemented as one or more intersection (collection)
of Slices.

Poly = {x : ∃ yᵢ, l ≤ a'x + b'yᵢ ≤ u for all slices, for all i ∈ secondary}
"""
#struct Poly
#    main::Set{Slice}
#end
Poly = Set{Slice}

"""
Creates closed Poly from Matrix-Vector representation.
"""
function Set{Slice}(A,l,u)
    @assert length(l) == length(u) == size(A,1)
    Poly(Slice(A[i,:], l[i], u[i]) for i in 1:length(l)) 
end
function Set{Slice}(A,l,u,rl,ru)
    @assert length(l) == length(u) == size(A,1) == length(rl) == length(ru)
    Poly(Slice(A[i,:], l[i], u[i],rl[i], ru[i]) for i in 1:length(l)) 
end

"""
Returns Matrix-vector Polyhedron form as vector
[A, l, u, rl, ru]
"""
function vectorize(p::Poly)
    # Warning: make sure that iterating over p gives same order each time
    A = reduce(vcat, (s.a' for s in p))
    l = reduce(vcat, ([s.l,] for s in p))
    u = reduce(vcat, ([s.u,] for s in p))
    rl = reduce(vcat, ([s.rl,] for s in p))
    ru = reduce(vcat, ([s.ru,] for s in p))
    (A,l,u,rl,ru)
end

function simplify(p::Poly; tol=1e-6)
    Keep = Dict{SparseVector{Float64,Int64}, Tuple{Float64,Float64,Relation,Relation}}()
    for s in p
        exists = false
        for (k,v) in Keep
            if isapprox(k, s.a; atol=tol)
                if v[1] > s.l + tol
                    l = v[1]
                    rl = v[3]
                elseif s.l > v[1] + tol
                    l = s.l
                    rl = s.rl
                else
                    l = 0.5*(v[1]+s.l)
                    rl = (v[3] == <) ? (<) : s.rl
                end
                if v[2] < s.u - tol
                    u = v[2]
                    ru = v[4]
                elseif s.u < v[2] - tol
                    u = s.u
                    ru = s.ru
                else
                    u = 0.5*(v[2]+s.u)
                    ru = (v[4] == <) ? (<) : s.ru
                end
                Keep[k] = (l,u,rl,ru)
                exists = true
                break
            end
        end
        if !exists && norm(s.a) > tol
            Keep[s.a] = (s.l, s.u, s.rl, s.ru)
        end
    end
    Poly(Slice(k, v...) for (k,v) in Keep)
end

"""
Utility method for easily generating random low-dimensional polyhedra.
"""
function Random.rand(rng::AbstractRNG, ::Random.SamplerType{Poly})
    m = rand(rng, 2:5)
    n = rand(rng, 2:5)
    A = sprandn(rng, n, m, 0.5)
    l = randn(rng, n)
    u = randn(rng, n)
    rl = rand(rng, Relation, n)
    ru = rand(rng, Relation, n)
    inconsistent = u .< l
    l[inconsistent] = u[inconsistent]
    rl[inconsistent] .= (≤)
    ru[inconsistent] .= (≤)
    Poly(A, l, u, rl, ru)
end

function random_polys_of_dim(N, m)
    polys = map(1:N) do i
        n = rand(2:4)
        A = sprandn(n, m, 0.75)
        l = randn(n)
        u = randn(n)
        rl = rand(Relation, n)
        ru = rand(Relation, n)
        inconsistent = u .< l
        l[inconsistent] = u[inconsistent] .- 2.0
        rl[inconsistent] .= (≤)
        ru[inconsistent] .= (≤)
        Poly(A, l, u, rl, ru)
    end
end

function embedded_dim(p::Poly)
    length(first(p).a) 
end

"""
Identify which constraints have non-infinite, open extents.
"""
function open_bounds(l, u, rl, ru)
    (; open_low=(rl .== (<)) .&& .!(isinf.(l)), open_hi=(ru .== (<)) .&& .!(isinf.(u)))
end

"""
Form the closure of the polyhedron.
"""
function closure(p::Poly)
    Poly(closure(s) for s in p)
end

"""
Convert to Polyhedra.jl polyhedron object. #TODO should probably use this format everywhere

NOTE assumes that p is a closed polyhedron. 
"""
function get_Polyhedron(p::Poly; tol=1e-6)
    hrep_poly = mapreduce(∩, p) do s
        cons = []
        if isapprox(s.l, s.u; atol=tol)
            push!(cons, Polyhedra.HyperPlane(s.a, s.u))
        else
            if !isinf(s.l)
                push!(cons, Polyhedra.HalfSpace(-s.a, -s.l))
            end
            if !isinf(s.u)
                push!(cons, Polyhedra.HalfSpace(s.a, s.u))
            end
        end
        reduce(∩, cons)
    end
    #Polyhedra.polyhedron(hrep_poly)
end

"""
Get vertices of poly.

This is probably not efficient for most polys of large size and low implicit dimension.
"""
function get_verts(p; tol=1e-6)
    hrep = get_Polyhedron(p; tol)
    vrep = Polyhedra.doubledescription(hrep)
    if length(vrep.points.points) == 0
        @infiltrate
    end
    #@assert length(vrep.rays.rays) == length(vrep.rays.lines.lines) == 0 # might not be necessary to check this
    (; V = vrep.points.points, R = vrep.rays.rays, L = vrep.rays.lines.lines)
end

"""
Project the poly into lower embedded dimension.
"""
function project(p::P, keep_dims; tol=1e-6) where P<:Union{Poly, Poly}
    hrep = get_Polyhedron(simplify(p); tol)
    vrep = Polyhedra.doubledescription(hrep)
    poly = Polyhedra.polyhedron(vrep)
    projected = Polyhedra.project(poly, keep_dims)
    local hrep
    try
        hrep = Polyhedra.doubledescription(projected.vrep)
    catch e 
        @infiltrate
    end
    AUi = mapreduce(vcat, Polyhedra.halfspaces(hrep); init=zeros(0, length(keep_dims)+1)) do hs
        [hs.a' hs.β]
    end
    AUe = mapreduce(vcat, Polyhedra.hyperplanes(hrep); init=zeros(0,length(keep_dims)+1)) do hp
        [hp.a' hp.β]
    end
    ni = size(AUi,1)
    ne = size(AUe,1)
    A = sparse([AUi[:,1:end-1]; AUe[:,1:end-1]])
    l = Vector{Float64}([fill(-Inf, ni); AUe[:,end]])
    u = Vector{Float64}([AUi[:,end]; AUe[:,end]])
    P(A,l,u) 
end

"""
Slice the poly given the specified values. 

Note that "slice" in this context is not to be confused with the custom type "Slice" which 
Poly objects are comprised of.
"""
function poly_slice(poly::Poly, x::Vector{Union{Float64, Missing}})
    keep_dims = ismissing.(x)
    slice_dims = .!keep_dims 
    sliced_poly = map(collect(poly)) do S
        slice_amount = S.a[slice_dims]'*x[slice_dims]
        a = S.a[keep_dims]
        l = S.l - slice_amount
        u = S.u - slice_amount
        Slice(a,l,u,S.rl,S.ru)
    end |> Set
end

"""
Return true if polyhedron is empty.
"""
function Base.isempty(poly::Poly; tol=1e-4, debug=false)
    m = OSQP.Model()
    n = length(poly)
    n == 0 && return false
    (A,l,u,rl,ru) = vectorize(poly)
    d = size(A,2)

    AA = [[A; A] zeros(2n); zeros(1,d+1)]
    (; open_low, open_hi) = open_bounds(l,u,rl,ru)
    l = [l; fill(-Inf, n); 0]
    u = [fill(Inf, n); u; 1]
    AA[1:n,end] = -1.0 * open_low
    AA[n+1:2n, end] = 1.0 * open_hi
    AA[end,end] = 1.0 
    
    OSQP.setup!(m;
                q=[zeros(d); -1.0],
                A=sparse(AA),
                l,
                u,
                verbose=false,
                polish=true,
                eps_abs=1e-8,
                eps_rel=1e-8)
    res = OSQP.solve!(m)
    @infiltrate debug
    if abs(res.info.status_val) == 3
        return true
    else
        t = -res.info.obj_val
        return t < tol
    end
end

"""
Identify constraints which have implicitly equivalent upper and lower bounds.
"""
function implicit_bounds(poly::Poly; tol=1e-4, debug=false)
    m = OSQP.Model()
    n = length(poly)
    (A,l,u,rl,ru) = vectorize(poly)
    A = sparse(A)
    implicit = fill(false, n)
    vals = fill(Inf, n)
    for i = n:-1:1
        if isapprox(l[i], u[i]; atol=tol)
            implicit[i] = true
            vals[i] = 0.5*(l[i] + u[i]) # trying to avoid tolerance issues
            continue
        else
            OSQP.setup!(m;
                        q = Vector{Float64}(A[i,:]),
                        A,
                        l,
                        u,
                        verbose=false,
                        polish=true,
                        eps_abs=1e-8,
                        eps_rel=1e-8)
            res = OSQP.solve!(m)
            if abs(res.info.status_val) == 3 # primal infeasible / empty set
                error("Empty set")
            end
            if abs(res.info.status_val) == 4 # dual infeasible / unbounded
                val_low = -Inf
            else
                val_low = res.info.obj_val
            end
            OSQP.setup!(m;
                        q = Vector{Float64}(-A[i,:]),
                        A,
                        l,
                        u,
                        verbose=false,
                        polish=true,
                        eps_abs=1e-8,
                        eps_rel=1e-8)
            res = OSQP.solve!(m)
            if abs(res.info.status_val) == 4 # dual infeasible / unbounded
                val_hi = Inf
            else
                val_hi = -res.info.obj_val
            end
            @infiltrate debug

            implicit[i] = isapprox(val_low, val_hi; atol=tol)
            implicit[i] && (vals[i] = 0.5*(val_hi+val_low))
        end
    end
    (; implicit, vals)
end

"""
Return the intrinsic dimension of the polyhedron.
"""
function intrinsic_dim(p::Poly; tol=1e-4, debug=false)
    local implicit, vals
    try
        (; implicit, vals) = implicit_bounds(p; tol, debug)
    catch e
        return 0 # TODO this is a hack... assuming that primal inf check only fails if intrinsic dim is 0... probably not true
    end
    (A,l,u,rl,ru) = vectorize(p)
    Aim = A[implicit,:]
    @infiltrate debug
    intrinsic_dim = embedded_dim(p) - rank(Aim)
end

"""
Return true if x is an element of p. Assumes closed poly.
"""
function Base.in(x::Vector{Float64}, p::Poly; tol=1e-6, debug=false)
    d = embedded_dim(p)
    n = length(x)
    @infiltrate debug
    if n == d
        return all( x in S for S in p )
    else
        (A,l,u,rl,ru) = vectorize(p)
        inds = collect(1:n)
        uninds = collect(n+1:d)
        Ap = A[:, inds]
        Ax = A[:,uninds]
        lx = length(uninds)
        m = OSQP.Model()
        OSQP.setup!(m; 
                    P=sparse(I,lx,lx), 
                    q = zeros(lx), 
                    A=sparse(Ax), 
                    l=l-Ap*x[inds], 
                    u=u-Ap*x[inds], 
                    polish=true, 
                    eps_abs=1e-8, 
                    eps_rel=1e-8, 
                    verbose=false)
        res = OSQP.solve!(m)
        sv = abs(res.info.status_val)
        @infiltrate (sv != 3 && sv != 1)
        return abs(res.info.status_val) != 3
    end
end
function Base.in(x::Vector{Float64}, s::Slice; tol=1e-6)
    ax = s.a'*x
    s.rl(s.l-tol, ax) && s.ru(ax-tol, s.u)
end

"""
Intersect p with polys.
"""
function poly_intersect(p::Poly, ps::Poly...)
    d = embedded_dim(p)
    @assert all(embedded_dim(psi) == d for psi in ps)
    union(p, ps...)
end

"""
Union of polyhedra.
"""
PolyUnion = Vector{Poly}

"""
Return true if x is an element of pu.
"""
function Base.in(x::Vector{Float64}, pu::PolyUnion)
    # TODO add tolerances
    any(x ∈ p for p in pu) 
end

"""
Complement of a polyhedra.
"""
function complement(s::Slice)
    out = PolyUnion()
    if !isinf(s.l)
        push!(out, Poly((Slice(s.a, -Inf, s.l, <, complement(s.rl)),)))
    end
    if !isinf(s.u)
        push!(out, Poly((Slice(s.a, s.u, Inf, complement(s.ru), <),)))
    end
    out
end
function complement(p::Poly)
    reduce(vcat, (complement(s) for s in p))
end

"""
Intersect unions of polyhedra.
"""
function poly_intersect(p::PolyUnion, ps::PolyUnion...)
    unions = (poly_intersect(subpieces...) for subpieces in Iterators.product(p, ps...))
end
