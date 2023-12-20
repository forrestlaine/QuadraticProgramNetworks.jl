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
returns (true iff x is lexico positive, magnitude of first non-zero value)
"""
function lexico_positive(x::Vector{Float64})
    return (first(x) ≥ 0, abs(first(x)))
end
function lexico_positive(x::SparseVector)
    (I, V) = findnz(x)
    t = argmin(I)
    return (V[t] ≥ 0, abs(V[t]))
end

function get_lexico_ordering(A)
    order = []
    for j = 1:size(A,2)
        for i = 1:size(A,1)
            a = A[i,:]
            (I, V) = findnz(a)
            if length(I) == 0 
                if j == 1
                    push!(order, i)
                end
                continue
            end
            if minimum(I) == j
                push!(order, i)
                continue
            end
        end
    end
    order
end


"""
A label which can be used to identify where halfspace constraints were introduced
in the QPN.
"""
struct HalfspaceLabel
    level::Int
    subpiece_index::Int
    comp_index::Int
    bound_index::Int
end

"""
A normalized slice.

S := { x : l ⋈ₗ a'x ⋈ᵤ u }

a ∈ ℝⁿ is required to be lexico-positive and have 1.0 as the leading non-zero value.

"""
struct Slice
    a::SparseVector{Float64, Int64}
    l::Float64
    u::Float64
    rl::Relation
    ru::Relation
    il::Set{HalfspaceLabel}
    iu::Set{HalfspaceLabel}
    Slice(a,l,u,rl::Relation,ru::Relation,il=Set{HalfspaceLabel}(),iu=Set{HalfspaceLabel}(); tol=1e-8) = begin
        droptol!(a, tol)
        n = norm(a)
        if (n ≤ tol)
            new(spzeros(length(a)), l, u, rl, ru, il, iu)
        else
            l_pos, n = lexico_positive(a)
            if l_pos
                new(sparse(a/n), l/n, u/n, rl, ru, il, iu)
            else
                new(sparse(-a/n), -u/n, -l/n, ru, rl, il, iu)
            end
        end
    end
    Slice(a,l,u) = Slice(a,l,u,≤,≤)
    Slice(a,l,u,il::Set{HalfspaceLabel},iu::Set{HalfspaceLabel}) = Slice(a,l,u,≤,≤,il,iu)
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
abstract type Poly end

struct BasicPoly <: Poly
    poly::Set{Slice}
end

struct ProjectedPoly <: Poly
    poly::BasicPoly
    parent::Poly
end

struct IntersectionPoly <: Poly
    polys::Vector{Poly}
end

struct LabeledPoly
    poly::Poly
    labels
end

function Base.isequal(S1::Poly, S2::Poly)
    isequal(Set(collect(S1)), Set(collect(S2)))
end
function Base.hash(S::Poly, h::UInt)
    hash(Set(collect(S)), h)
end

"""
Creates closed Poly from Matrix-Vector representation.
"""
function Poly(A,l,u)
    @assert length(l) == length(u) == size(A,1)
    BasicPoly(Set(Slice(A[i,:], l[i], u[i]) for i in 1:length(l)))
end
function Poly(A,l,u,rl::Vector{Relation},ru::Vector{Relation})
    @assert length(l) == length(u) == size(A,1) == length(rl) == length(ru)
    BasicPoly(Set(Slice(A[i,:], l[i], u[i], rl[i], ru[i]) for i in 1:length(l)))
end
function Poly(A,l,u,il::Vector{Set{HalfspaceLabel}},iu::Vector{Set{HalfspaceLabel}})
    @assert length(l) == length(u) == size(A,1) == length(il) == length(iu)
    BasicPoly(Set(Slice(A[i,:], l[i], u[i], il[i], iu[i]) for i in 1:length(l)))
end
function Poly(S)
    BasicPoly(Set(S))
end

function Base.eltype(poly::Poly)
    Slice
end
function Base.IteratorSize(poly::Poly)
    Base.HasLength()
end

function Base.length(poly::Union{BasicPoly,ProjectedPoly})
    length(poly.poly)
end
function Base.length(poly::IntersectionPoly)
    sum(length(p) for p in poly.polys)
end
function Base.iterate(poly::Union{BasicPoly,ProjectedPoly})
    iterate(poly.poly)
end
function Base.iterate(poly::Union{BasicPoly,ProjectedPoly}, state)
    iterate(poly.poly, state)
end
function Base.iterate(poly::IntersectionPoly)
    if length(poly.polys) == 0
        return nothing
    else
        iterate(poly, (; current=1, inner_state=nothing) )
    end
end
function Base.iterate(poly::IntersectionPoly, state)
    if state.current > length(poly.polys)
        return nothing
    end
    if isnothing(state.inner_state)
        ret = iterate(poly.polys[state.current])
    else
        ret = iterate(poly.polys[state.current], state.inner_state)
    end
    if isnothing(ret)
        return iterate(poly, (; current=state.current+1, inner_state=nothing))
    else
        return (ret[1], (; current=state.current, inner_state=ret[2]))
    end
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
    (;A,l,u,rl,ru)
end

function has_parent(p::BasicPoly, i)
    false
end
function has_parent(p::ProjectedPoly, i)
    true
end
function has_parent(p::IntersectionPoly, i)
    id = 0
    for pp in p.polys
        len = length(pp)
        if id+1 ≤ i ≤ id+len
            return has_parent(pp, i-id)
        else
            id += len
        end
    end
end
function get_parent(p::ProjectedPoly, i)
    p.parent
end
function get_parent(p::IntersectionPoly, i)
    id = 0
    for pp in p.polys
        len = length(pp)
        if id+1 ≤ i ≤ id+len
            return pp.parent
        else
            id += len
        end
    end
end

function simplify(p::BasicPoly; tol=1e-6)
    #try
    #    #local implicilty_equality, vals
    #    (implicitly_equality, vals) = implicit_bounds(p; tol, debug)
    #    (;A,l,u,rl,ru,il,iu) = vectorize(p)
    #    l[implicitly_equality] = u[implicitly_equality] = vals
    #    p = Poly(A,l,u,rl,ru,il,iu)
    #catch e
    #end
    Keep = Dict{SparseVector{Float64,Int64}, Tuple{Float64,Float64,Relation,Relation,Set{HalfspaceLabel}, Set{HalfspaceLabel}}}()
    for s in p
        exists = false
        for (k,v) in Keep
            if isapprox(k, s.a; atol=tol)
                if v[1] > s.l + tol
                    l = v[1]
                    rl = v[3]
                    il = v[5]
                elseif s.l > v[1] + tol
                    l = s.l
                    rl = s.rl
                    il = s.il
                else
                    l = 0.5*(v[1]+s.l)
                    rl = (v[3] == <) ? (<) : s.rl
                    il = s.il ∪ v[5]
                end
                if v[2] < s.u - tol
                    u = v[2]
                    ru = v[4]
                    iu = v[6]
                elseif s.u < v[2] - tol
                    u = s.u
                    ru = s.ru
                    iu = s.iu
                else
                    u = 0.5*(v[2]+s.u)
                    ru = (v[4] == <) ? (<) : s.ru
                    iu = s.iu ∪ v[6]
                end
                Keep[k] = (l,u,rl,ru,il,iu)
                exists = true
                break
            end
        end
        if !exists && norm(s.a) > tol
            Keep[s.a] = (s.l, s.u, s.rl, s.ru, s.il, s.iu)
        end
    end
    Poly(Slice(k, v...) for (k,v) in Keep)
end
function simplify(p::ProjectedPoly; tol=1e-6)
    ProjectedPoly(simplify(p.poly; tol), p.parent)
end
function simplify(p::IntersectionPoly; tol=1e-6)
    IntersectionPoly([simplify(poly; tol) for poly in p.polys])
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
#function closure(p::Poly)
#    Poly(closure(s) for s in p)
#end
function closure(p::BasicPoly)
    BasicPoly(Set(closure(s) for s in p))
end
function closure(p::ProjectedPoly)
    ProjectedPoly(closure(p.poly), p.parent)
end
function closure(p::IntersectionPoly)
    IntersectionPoly([closure(sub_p) for sub_p in p.polys])
end

"""
Determine whether P1 ⊆ P2, i.e. if P1 is a subset of P2.
"""
function Base.issubset(P1::Poly, P2::Poly; tol=1e-6)
    (; A, l, u) = vectorize(P1)
    A1 = A; l1 = l; u1 = u;
    (; A, l, u) = vectorize(P2)
    A2 = A; l2 = l; u2 = u;

    m = length(l2)
    for i = 1:m
        for (bound, dir) in zip((l2[i], u2[i]), (1.0, -1.0))
            if isfinite(bound)
                model = OSQP.Model()
                OSQP.setup!(model; 
                            q = dir*collect(A2[i,:]), 
                            A=sparse(A1), 
                            l=l1, 
                            u=u1, 
                            verbose=false, 
                            polish=true, 
                            eps_abs=tol, 
                            eps_rel=tol)
                ret = OSQP.solve!(model)
                if ret.info.status_val != 1
                    return false # appears unbounded below
                elseif ret.info.obj_val < dir*bound - tol
                    return false
                end
            end
        end
    end
    return true
end


"""
Convert to Polyhedra.jl polyhedron object. #TODO should probably use this format everywhere

NOTE assumes that p is a closed polyhedron. 
"""
function get_Polyhedron_hrep(p::Poly; tol=1e-6)

    hyperplanes = Polyhedra.HyperPlane{Float64, SparseVector{Float64, Int64}}[]
    halfspaces = Polyhedra.HalfSpace{Float64, SparseVector{Float64, Int64}}[]
    for s in p
        if isapprox(s.l, s.u; atol=tol) 
            push!(hyperplanes, Polyhedra.HyperPlane(s.a, s.u))
        else
            if !isinf(s.l)
                push!(halfspaces, Polyhedra.HalfSpace(-s.a, -s.l))
            end
            if !isinf(s.u)
                push!(halfspaces, Polyhedra.HalfSpace(s.a, s.u))
            end
        end
    end
    hr = Polyhedra.hrep(hyperplanes, halfspaces)
end

"""
Get vertices of poly.

This is probably not efficient for most polys of large size and low implicit dimension.
"""
function get_verts(p; tol=1e-6)
    hr = get_Polyhedron_hrep(p; tol)
    vr = Polyhedra.doubledescription(hr)
    if length(vr.points.points) == 0
        error("There should be at least one vertex in a non-empty pointed cone")
    end
    (; V = vr.points.points, R = vr.rays.rays, L = vr.rays.lines.lines)
end

"""
Converts a Polyhedra.VRep object to a QPN Poly object
"""
function vrep_to_poly(vr)
    polyh = polyhedron(vr)
    #proj_hr = hrep(projected)
    HS = Polyhedra.halfspaces(polyh)
    HP = Polyhedra.hyperplanes(polyh)
    dim = 0
    if length(HS) == 0
        AUi = nothing
        ni = 0
    else
        AUi = mapreduce(vcat, HS) do hs
            [hs.a' hs.β]
        end
        ni = size(AUi,1)
        dim = size(AUi, 2)
    end
    if length(HP) == 0
        AUe = nothing
        ne = 0
    else
        AUe = mapreduce(vcat, HP) do hp
            [hp.a' hp.β]
        end
        ne = size(AUe,1)
        dim = size(AUe, 2)
    end
    if dim == 0 
        error("Vrep appears to be empty")
    else
        if isnothing(AUi)
            AUi = zeros(0, dim)
        end
        if isnothing(AUe)
            AUe = zeros(0, dim)
        end
    end
    A = sparse([AUi[:,1:end-1]; AUe[:,1:end-1]])
    l = Vector{Float64}([fill(-Inf, ni); AUe[:,end]])
    u = Vector{Float64}([AUi[:,end]; AUe[:,end]])
    Poly(A,l,u)
end

"""
Project the poly into lower embedded dimension.
"""
function project(p::Poly, keep_dims; tol=1e-6)
    hr = get_Polyhedron_hrep(p; tol)
    poly = Polyhedra.polyhedron(hr)
    vr = vrep(poly)

    N = embedded_dim(p)

    Pmat = sparse(I, N, N)
    Pmat = Pmat[keep_dims, :]

    projected_points::Vector{Vector{Float64}} = map(point->Pmat*point, Polyhedra.points(vr))
    projected_rays::Vector{Polyhedra.Ray{Float64, Vector{Float64}}} = map(ray->Pmat*ray, Polyhedra.rays(vr))
    projected_lines::Vector{Polyhedra.Line{Float64, Vector{Float64}}} = map(line->Pmat*line, Polyhedra.lines(vr))

    proj_vr = vrep(projected_points, projected_lines, projected_rays)
    projected = vrep_to_poly(proj_vr)
    ProjectedPoly(projected, p)
end


"""
Slice the poly given the specified values. 

Note that "slice" in this context is not to be confused with the custom type "Slice" which 
Poly objects are comprised of.
"""
function poly_slice(poly::BasicPoly, x::Vector{Union{Float64, Missing}})
    keep_dims = ismissing.(x)
    slice_dims = .!keep_dims 
    sliced_poly = map(collect(poly)) do S
        slice_amount = S.a[slice_dims]'*x[slice_dims]
        a = S.a[keep_dims]
        l = S.l - slice_amount
        u = S.u - slice_amount
        Slice(a,l,u,S.rl,S.ru,S.il,S.iu)
    end |> Set |> BasicPoly
end
function poly_slice(poly::ProjectedPoly, x::Vector{Union{Float64, Missing}})
    ProjectedPoly(poly_slice(poly.poly, x), poly.parent)
end
function poly_slice(poly::IntersectionPoly, x::Vector{Union{Float64, Missing}})
    IntersectionPoly([poly_slice(p, x) for p in poly.polys])
end

function exemplar(poly::Poly; tol=1e-4, debug=false)
    m = OSQP.Model()
    n = length(poly)
    n == 0 && return (; empty=false, example=nothing)
    (; A,l,u,rl,ru) = vectorize(poly)
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
    t = -res.info.obj_val
    empty = abs(res.info.status_val) == 3 || t < tol
    example = empty ? nothing : res.x[1:end-1]
    (; empty, example)
end

"""
Return true if polyhedron is empty.
"""
function Base.isempty(poly::Poly; tol=1e-4, debug=false)
    (; empty) = exemplar(poly; tol, debug)
    return empty
end

"""
Identify constraints which have implicitly equivalent upper and lower bounds.
"""
function implicit_bounds(poly::Poly; tol=1e-4, debug=false)
    m = OSQP.Model()
    n = length(poly)
    (; A,l,u,rl,ru) = vectorize(poly)
    A = sparse(A)
    implicitly_equality = fill(false, n)
    vals = fill(Inf, n)
    for i = n:-1:1
        if isapprox(l[i], u[i]; atol=tol)
            implicitly_equality[i] = true
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

            implicitly_equality[i] = isapprox(val_low, val_hi; atol=tol)
            implicitly_equality[i] && (vals[i] = 0.5*(val_hi+val_low))
        end
    end
    (; implicitly_equality, vals)
end

"""
Return the intrinsic dimension of the polyhedron.
"""
function intrinsic_dim(p::Poly; tol=1e-4, debug=false)
    local implicitly_equality, vals
    try
        (; implicitly_equality, vals) = implicit_bounds(p; tol, debug)
    catch e
        return 0 # TODO this is a hack... assuming that primal inf check only fails if intrinsic dim is 0... probably not the case
    end
    (; A,l,u,rl,ru) = vectorize(p)
    Aim = A[implicitly_equality,:]
    @infiltrate debug
    intrinsic_dim = embedded_dim(p) - rank(Aim)
end

function eliminate_variables(p::Poly, indices, xz; debug=false)

    # TODO WARNING halfspace labels are deleted when using this method, need to
    # fix

    elim_inds = indices
    keep_inds = setdiff(1:embedded_dim(p), elim_inds)
    if keep_inds == 1:embedded_dim(p)
        return p
    end
    local implicitly_equality, vals
    try
        (; implicitly_equality, vals) = implicit_bounds(p)
    catch e
        @error e
        return p
    end
    (; A,l,u,rl,ru) = vectorize(p)

    inequality = .!implicitly_equality
    Ae_elim = A[implicitly_equality,elim_inds]
    Ae_keep = A[implicitly_equality,keep_inds]
    Ai_elim = A[inequality,elim_inds]
    Ai_keep = A[inequality,keep_inds]

    xz_keep = xz[keep_inds]
    xz_elim = xz[elim_inds]

    rhs = vals[implicitly_equality]

    F = qr(Ae_elim)

    if rank(F) < size(Ae_elim, 2)
        @warn "Not enough constraints to eliminate. rank is $(rank(Ae_elim)) and size is $(size(Ae_elim))"

        qr_keep = []
        offset = 0
        for i = 1:size(Ae_elim, 2)
            if i+offset > size(Ae_elim, 2)
                break
            end
            if !iszero(F.R[i,i+offset])
                push!(qr_keep, i)
            else
                offset += 1
            end
        end
        new_elim_inds = elim_inds[F.pcol[qr_keep]]
        new_keep_inds = [keep_inds; setdiff(elim_inds, new_elim_inds)]

        @info "Can't eliminate following inds: $(setdiff(elim_inds, new_elim_inds))"
        
        keep_inds = new_keep_inds
        elim_inds = new_elim_inds

        Ae_elim = A[implicitly_equality,elim_inds]
        Ae_keep = A[implicitly_equality,keep_inds]
        Ai_elim = A[inequality,elim_inds]
        Ai_keep = A[inequality,keep_inds]

        xz_keep = xz[keep_inds]
        xz_elim = xz[elim_inds]

        @assert rank(Ae_elim) == size(Ae_elim, 2)

    end
    # x2 = (Ae_elim)† * (rhs - Ae_keep * x1)
    Ad = (Ae_elim'*Ae_elim)\Matrix(Ae_elim') |> sparse
    P = (I - Ae_elim*Ad)
    Ae = P*Ae_keep
    be = P*rhs
    Ai = Ai_keep - Ai_elim*Ad*Ae_keep
    ci = Ai_elim*Ad*rhs
    ui = u[inequality] - ci
    li = l[inequality] - ci

    @infiltrate debug

    return Poly([Ae;Ai], 
                [be;li], 
                [be;ui], 
                [rl[implicitly_equality]; rl[inequality]],
                [ru[implicitly_equality]; ru[inequality]])
end
    

"""
Return true if x is an element of p. Assumes closed poly.
"""
function Base.in(x::Vector{Float64}, p::Poly; tol=1e-6, debug=false)
    d = embedded_dim(p)
    n = length(x)
    @infiltrate debug
    if n == d
        return all( in(x, S; tol) for S in p )
    else
        (; A,l,u,rl,ru) = vectorize(p)
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
                    eps_abs=tol, 
                    eps_rel=tol, 
                    verbose=false)
        res = OSQP.solve!(m)
        sv = abs(res.info.status_val)
        @infiltrate (sv ∉ (1,2,3))
        return abs(res.info.status_val) != 3
    end
end
function Base.in(x::Vector{Float64}, s::Slice; tol=1e-6)
    ax = s.a'*x
    s.rl(s.l-tol, ax) && s.ru(ax-tol, s.u)
end

"""
Union of polyhedra.
"""
struct PolyUnion
    polys::Vector{Poly}
end
function Base.eltype(pu::PolyUnion)
    Poly
end
function Base.IteratorSize(pu::PolyUnion)
    Base.HasLength()
end
function Base.length(pu::PolyUnion)
    length(pu.polys)
end
function Base.iterate(pu::PolyUnion)
    iterate(pu.polys)
end
function Base.iterate(pu::PolyUnion, state)
    iterate(pu.polys, state)
end
function Base.vcat(pus::PolyUnion...)
    PolyUnion(reduce(vcat, pu.polys for pu in pus))
end
function Base.getindex(pu::PolyUnion, i)
    pu.polys[i]
end
function Base.firstindex(pu::PolyUnion)
    firstindex(pu.polys)
end
function Base.lastindex(pu::PolyUnion)
    lastindex(pu.polys)
end

function remove_subsets(pu::PolyUnion)
    is_subset = zeros(Bool, length(pu))
    Threads.@threads for i in 1:length(pu)
        if any(i ≠ j && !is_subset[j] && pu[i] ⊆ p for (j,p) in enumerate(pu))
            is_subset[i] = true
        end
    end
    proper_sets = pu[.!is_subset]
    sum(is_subset) > 0 && @info "Removed $(sum(is_subset)) redundant sets."
    return PolyUnion(proper_sets)
end
function remove_subsets(pu::Nothing)
    nothing
end

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
    out = BasicPoly[]
    if !isinf(s.l)
        push!(out, Poly((Slice(s.a, -Inf, s.l, <, complement(s.rl)),)))
    end
    if !isinf(s.u)
        push!(out, Poly((Slice(s.a, s.u, Inf, complement(s.ru), <),)))
    end
    PolyUnion(out)
end
function complement(p::Poly)
    reduce(vcat, (complement(s) for s in p))
end

"""
Intersect p with polys.
"""

function poly_intersect(p::Poly, ps::Poly...)
    if p isa IntersectionPoly
        polys = p.polys
    else
        polys = [p]
    end
    for pp in ps
        if pp isa IntersectionPoly
            append!(polys, pp.polys)
        else
            push!(polys, p)
        end
    end
    IntersectionPoly(polys)
end
function poly_intersect(p::Union{BasicPoly, ProjectedPoly}...)
    d = embedded_dim(first(p))
    @assert all(embedded_dim(psi) == d for psi in p)
    IntersectionPoly(vcat(p...))
end

function poly_intersect(p::Union{BasicPoly, ProjectedPoly}, ip::IntersectionPoly)
    @assert embedded_dim(p) == embedded_dim(ip)
    IntersectionPoly([p; ip.polys])
end

function poly_intersect(p::IntersectionPoly...)
    @assert allequal(embedded_dim(pp) for pp in p)
    all_polys = mapreduce(vcat, p) do pp
        pp.polys
    end
    IntersectionPoly(all_polys)
end

"""
Intersect unions of polyhedra. (returns an iterator)
"""
function poly_intersect(p::PolyUnion, ps::PolyUnion...)
    unions = (poly_intersect(subpieces...) for subpieces in Iterators.product(p, ps...))
end

function convex_hull2(pu::PolyUnion)
    VV = Set{QuantizedVector}() 
    LL = Set{QuantizedVector}() 
    RR = Set{QuantizedVector}()
    for p in pu
        (; V, L, R) = get_verts(p)
        foreach(v->push!(VV,QuantizedVector(v=v)), V)
        foreach(l->push!(LL,QuantizedVector(v=l.a)), L)
        foreach(r->push!(RR,QuantizedVector(v=r.a)), R)
    end
    VV = map(collect(VV)) do v
        v.v
    end
    LL = map(collect(LL)) do l
        Polyhedra.Line(l.v)
    end
    RR = map(collect(RR)) do r
        Polyhedra.Ray(r.v)
    end
    vr = vrep(VV,LL,RR)
    vrep_to_poly(vr)
end
function convex_hull(pu::PolyUnion; tol=1e-6)
    polyhedra = map(pu) do p 
        hr = get_Polyhedron_hrep(simplify(p); tol)
        poly = Polyhedra.polyhedron(hr)
        vr = vrep(poly)
        poly
    end
    hull = Polyhedra.convexhull(polyhedra...)
    Polyhedra.removevredundancy!(hull)
    vr = vrep(hull)
    vrep_to_poly(vr)
end

"""
Determine whether P1 ⊆ P2, i.e. if P1 is a subset of P2.
"""
function Base.issubset(P1::Poly, P2::PolyUnion; tol=1e-6)
    @warn "Determining subset relations are hard for unions of polyhedra. This method is therefore not guaranteed to be correct. Return values of 'true' are correct, whereas a 'false' return may not actually imply a negative result."
    any(P1 ⊆ P for P in P2)
end
