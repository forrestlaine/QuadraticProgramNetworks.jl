
"""
Not-necessarily closed polyhedron:
{ x : lᵢ ⋈ᵢ (aᵢ'x) ⋈ᵢ uᵢ, i ∈ 1...m }
⋈ᵢ ∈ { ≤, < }
(lᵢ≤uᵢ) ∈ ℝ̄ (extended reals)
"""
Base.@kwdef struct Poly
    A::SparseMatrixCSC{Float64, Int64}
    l::Vector{Float64}
    u::Vector{Float64}
    rl::Vector{Relation}
    ru::Vector{Relation}
    Poly(A,l,u,rl,ru) = begin
        (size(A,1) == length(l) == length(u) == 
         length(rl) == length(ru)) ? 
        new(A,l,u,rl,ru) : error("Inconsistent sizes.")
    end
end

function Poly(p::Poly2)
    Poly(vectorize(p)...)
end

"""
ℝᵈ
"""
function Poly(d::Int)
    Poly(spzeros(0, d), zeros(0), zeros(0), fill(≤, 0), fill(≤, 0))
end

"""
Default constructor assumes closed polyhedron.
"""
function Poly(A::SparseMatrixCSC{Float64, Int64}, l::Vector{Float64}, u::Vector{Float64})
    cl = fill(≤, length(l))
    Poly(A, l, u, cl, cl)
end

"""
DON'T ACTUALLY DO THIS -- only for debugging convenience
"""
function Base.isequal(p1::Poly, p2::Poly) 
    p1.A == p2.A && p1.l == p2.l && p1.u == p2.u && p1.rl == p2.rl && p1.ru == p2.ru
end

"""
SEE ABOVE
"""
function Base.hash(p::Poly, h::UInt)
    hash(Set((p.A[i,:], p.l[i], p.u[i], p.rl[i], p.ru[i]) for i in 1:length(p)), h)
    #hash((p.A, p.l, p.u, p.rl, p.ru), h)
end
"""
Number of constraints comprising polyhedron.
"""
function Base.length(p::Poly)
    size(p.A,1)
end
#function Base.length(p::Poly2)
#    ALREADY IMPLEMENTED
#end

"""
The dimension of the space that the polyhedron is embedded in.
"""
function embedded_dim(p::Poly)
    size(p.A,2)
end

function open_bounds(p::Poly)
    (; open_low=(p.rl .== (<)) .&& .!(isinf.(p.l)), open_hi=(p.ru .== (<)) .&& .!(isinf.(p.u)))
end
function closure(p::Poly)
    cl = fill(≤, length(p))
    Poly(copy(p.A), copy(p.l), copy(p.u), cl, cl)
end
function get_Polyhedron(p::Poly; tol=1e-6)
    hrep_poly = mapreduce(∩, 1:length(p)) do i
        cons = []
        if isapprox(p.l[i], p.u[i]; atol=tol)
            push!(cons, Polyhedra.HyperPlane(p.A[i,:], p.l[i]))
        else
            if !isinf(p.l[i])
                push!(cons, Polyhedra.HalfSpace(-p.A[i,:], -p.l[i]))
            end
            if !isinf(p.u[i])
                push!(cons, Polyhedra.HalfSpace(p.A[i,:], p.u[i]))
            end
        end
        reduce(∩, cons)
    end
    poly = Polyhedra.polyhedron(hrep_poly)
    #Polyhedra.removehredundancy!(poly)
end
function inconsistent_poly(p)
    for i in 1:length(p)
        if p.l[i] > p.u[i]
            return true
        end
        if all(iszero.(p.A[i,:])) && (p.l[i] > 0 || p.u[i] < 0)
            return true
        end
    end
    false
end
function Base.isempty(poly::Poly; tol=1e-4, debug=false)
    inconsistent_poly(poly) && return true
    m = OSQP.Model()
    n = length(poly)
    d = embedded_dim(poly)


    AA = [[poly.A; poly.A] zeros(2n); zeros(1,d+1)]
    (; open_low, open_hi) = open_bounds(poly)
    l = [poly.l; fill(-Inf, n); 0]
    u = [fill(Inf, n); poly.u; 1]
    
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
function implicit_bounds(poly::Poly; tol=1e-4, debug=false)
    m = OSQP.Model()
    n = length(poly)
    A = sparse(poly.A)
    implicit = fill(false, n)
    vals = fill(Inf, n)
    for i = n:-1:1
        if isapprox(poly.l[i], poly.u[i]; atol=tol)
            implicit[i] = true
            vals[i] = 0.5*(poly.l[i] + poly.u[i]) # trying to avoid tolerance issues
            continue
        else
            OSQP.setup!(m;
                        q = Vector{Float64}(poly.A[i,:]),
                        A,
                        l=poly.l,
                        u=poly.u,
                        verbose=false,
                        polish=true,
                        eps_abs=1e-8,
                        eps_rel=1e-8)
            res = OSQP.solve!(m)
            if abs(res.info.status_val) == 3 # primal infeasible / empty set
                @infiltrate
                error("Empty set")
            end
            if abs(res.info.status_val) == 4 # dual infeasible / unbounded
                val_low = -Inf
            else
                val_low = res.info.obj_val
            end
            OSQP.setup!(m;
                        q = Vector{Float64}(-poly.A[i,:]),
                        A,
                        l=poly.l,
                        u=poly.u,
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
Make implicit constraint bounds explicit.
"""
function simplify!(p::Poly; tol=1e-4)
    (; implicit, vals) = implicit_bounds(p; tol)
    simplify!(p, implicit, vals)
end

"""
Make implicit constraint bounds explicit.
"""
function simplify!(p::Poly, implicit, vals)
    p.l[implicit] = vals[implicit]
    p.u[implicit] = vals[implicit]
    p.rl[implicit] .= ≤
    p.ru[implicit] .= ≤
end
function intrinsic_dim(p::Poly; tol=1e-4, debug=false)
    local implicit, vals
    try
        (; implicit, vals) = implicit_bounds(p; tol, debug)
    catch e
        @infiltrate
    end
    Aim = p.A[implicit,:]
    @infiltrate debug
    intrinsic_dim = embedded_dim(p) - rank(Aim)
end
function Base.in(x::Vector{Float64}, p::Poly; tol=1e-4)
    ax = p.A*x
    return all(p.rl[i](p.l[i]-tol, ax[i]) && p.ru[i](ax[i]-tol, p.u[i]) for i in 1:length(p))
end
function Base.intersect(p::Poly, polys::Poly...)
    A = vcat(p.A, (poly.A for poly in polys)...)
    l = vcat(p.l, (poly.l for poly in polys)...)
    u = vcat(p.u, (poly.u for poly in polys)...)
    rl = vcat(p.rl, (poly.rl for poly in polys)...)
    ru = vcat(p.ru, (poly.ru for poly in polys)...)
    Poly(A, l, u, rl, ru)
end

PolyUnion = Vector{Poly}
struct PolyUnionIterator
    len::Int
    polys
end

struct IntersectionTriplet
    primal::Bool
    ind1::Int
    ind2::Int
end
mutable struct PolyUnionIntersectionIterator
    L1::Vector{PolyUnionIterator}
    L2::Vector{PolyUnionIterator}
    state::Set{Set{IntersectionTriplet}}
end

function Base.length(pu::PolyUnionIterator)
    pu.len
end
function Base.length(pui::PolyUnionIntersectionIterator)
    prod(length(l1)+length(l2) for (l1,l2) ∈ zip(L1,L2))
end

Base.eltype(::Type{PolyUnionIterator}) = Poly2
Base.eltype(::Type{PolyUnionIntersectionIterator}) = Poly2

Base.iterate(pu::PolyUnionIterator) = iterate(pu.polys)
Base.iterate(pu::PolyUnionIterator, state) = iterate(pu.polys, state)
function complement(poly::Poly)
    n = length(poly)
    d = embedded_dim(poly)
    polys = []
    for i in 1:n
        if !isinf(poly.l[i])
            push!(polys, Poly(A=poly.A[i:i,:], l=[-Inf,], u=poly.l[i:i], rl = [<,], ru = complement.(poly.rl[i:i])))
        end
        if !isinf(poly.u[i])
            push!(polys, Poly(A=poly.A[i:i,:], l=poly.u[i:i], u=[Inf,], rl = complement.(poly.ru[i:i]), ru=[<,]))
        end
    end
    PolyUnion(polys)
end
function Base.show(io::IO, ::MIME"text/plain", poly::Poly)
    indent=get(io, :indent, 0)
    space = " "^indent
    n = length(poly)
    d = embedded_dim(poly)
    half = Int(ceil(n/2))
    println(io, space*"Polyhedron in ℝ^", d, " with ", n, " constraints.")
    for i = 1:n
        str = space*"%5.2f %2s | "
        args = [poly.l[i], poly.rl[i]]
        for j in 1:d
            if iszero(poly.A[i,j])
                str *= "  ⋅   "
            else
                str *= "%5.2f "
                push!(args, poly.A[i,j])
            end
        end
        if i == half
            str *= "| x %2s %5.2f"
        else
            str *= "|   %2s %5.2f"
        end
        push!(args, poly.ru[i])
        push!(args, poly.u[i])
        str *= "\n"
        format = Printf.Format(str)
        Printf.format(io, format, args...)
    end
end
function combine(regions, solutions, level_dim)
    all_combos = Iterators.product(zip(map(complement, regions), solutions)...)
    (_, all_combos) = Iterators.peel(all_combos) # skip first combo, corresponding to all region complements
    pieces = imap(all_combos) do combo_of_sets
        poly_intersect((combo_of_sets .|> collect)...) 
    end |> Iterators.flatten
    pieces = distinct(Iterators.filter(x->!isempty(x) && intrinsic_dim(x) ≥ level_dim, pieces))
end
#"""
#Given z solving avi(w)
#create the union of polyhedral solution regions 
#local to (z,w).
#"""
#function local_pieces(avi, z, w; debug=false, expansion=true, tol=1e-5)
#    n = length(z)
#    m = length(w)
#    J = comp_indices(avi, z, w)
#    
#    J2 = Set(J[2])
#    J4 = Set(J[4])
#    Ks = map(enumerate(Iterators.product(powerset(J[2]), powerset(J[4])))) do (e,(S2, S4))
#        C2 = setdiff(J2, Set(S2)) |> collect
#        C4 = setdiff(J4, Set(S4)) |> collect
#        Dict(1=>Set([J[1];C2]), 2=>Set([J[3];S2;S4]),3=>Set([J[5];C4]),4=>Set(J[6]))
#    end |> Set
#
#    extended_Ks = Ks
#    vertices = Vector{Float64}[]
#    processed_pieces = Poly[]
#    depth = 1 
#    while true
#        new_pieces = Iterators.filter(!isempty, (local_piece(avi,n,m,K) for K in extended_Ks)) |> collect
#        processed_pieces = [processed_pieces; new_pieces]
#        !expansion && break
#        new_verts = Vector{Float64}[]
#        for piece in new_pieces
#            local_verts = [[v;w] for v in get_verts(poly_slice(piece, [fill(missing, n); w])).V]
#            piece_verts = [v for v in local_verts if !(any(vv->isapprox(v, vv; atol=tol), [vertices; new_verts]))]
#            new_verts = [new_verts; piece_verts] 
#        end
#
#        extended_Ks = mapreduce(union, new_verts; init=[]) do v
#            JJ = comp_indices(avi, v[1:length(z)], w) 
#            J2 = Set(JJ[2])
#            J4 = Set(JJ[4]) 
#            Kv = map(enumerate(Iterators.product(powerset(JJ[2]), powerset(JJ[4])))) do (e,(S2, S4))
#                C2 = setdiff(J2, Set(S2)) |> collect
#                C4 = setdiff(J4, Set(S4)) |> collect
#                K = Dict(1=>Set([JJ[1];C2]), 2=>Set([JJ[3];S2;S4]),3=>Set([JJ[5];C4]),4=>Set(JJ[6]))
#            end |> Set
#        end |> (extended_Ks -> setdiff(extended_Ks, Ks))
#
#        if isempty(extended_Ks)
#            break
#        else
#            union!(Ks, extended_Ks)
#            vertices = [vertices; new_verts]
#            depth += 1
#        end
#    end
#    debug && println("Just generated ", length(processed_pieces), " solutions to AVI, exploring depth of ", depth, ".")
#    return processed_pieces
#end
#
#"""
#NEW VERSION:
#    Should find all polyhedral regions Rᵢ(w) s.t. z ∈ Rᵢ(w) ⟹  z solves AVI(w).
#    Regions should only be local, meaning for given solution pair (w,z), 
#    (z) ∈ Rᵢ(w)
#
#    P = { d : dᵢ = 0,             i ∈ J₁ ∪ J₅ ∪ J₆,
#              lᵢ ≤ zᵢ + dᵢ ≤ uᵢ,  i ∈ J₂ ∪ J₃ ∪ J₄ }
#    R = []
#    for (i, (S₂, S₄)) ∈ enumerate(product(powerset(J₂), powerset(J₄)))
#        C2 = setdiff(J2, Set(S2))
#        C4 = setdiff(J4, Set(S4))
#        Rᵢ = P ∩ { d : (Md)ⱼ ≥ 0, dⱼ = 0, j ∈ S₂,
#                       (Md)ⱼ = 0, dⱼ ≥ 0, j ∈ C₂, 
#                       (Md)ⱼ ≤ 0, dⱼ = 0, j ∈ S₄,
#                       (Md)ⱼ = 0, dⱼ ≤ 0, j ∈ C₄ }
#        push!(R, Rᵢ)
#    end
#
#OLD VERSION:
#J[1]    : Md ≥ 0 (I think this is wrong... should be free?)
#J[2]    : Md ≥ 0 
#J[3]    : Md = 0
#J[4]    : Md ≤ 0
#J[5]    : Md ≤ 0 (I think this is wrong... should be free?)
#all     : l ≤ z+d ≤ u
#        : (M+M')d = 0
#
#
#I think this is wrong...
#
## TODO probably inefficient initialization with Poly
#"""
#function single_piece(avi, z, w, J)
#    return nothing # depricated
#    A = [avi.M[[J[1];J[2]],:]; avi.M[J[3],:]; avi.M[[J[4]; J[5]],:]; avi.M+avi.M'; I(length(z))]
#    l = [zeros(length(J[1])+length(J[2])+length(J[3])); fill(-Inf, length(J[4])+length(J[5])); zeros(length(z)); avi.l-z]
#    u = [fill(Inf, length(J[1])+length(J[2])); zeros(length(J[3])+length(J[4])+length(J[5])); zeros(length(z)); avi.u-z]
#
#    meaningful = find_non_trivial(A,l,u)
#    P1 = Poly([A[meaningful,:] spzeros(sum(meaningful), length(w))], l[meaningful], u[meaningful])
#    
#    A = [avi.M[J[2],:]; avi.M[J[3],:]; avi.M[J[4],:]; avi.M+avi.M'; I(length(z))]
#    l = [zeros(length(J[2])+length(J[3])); fill(-Inf, length(J[4])); zeros(length(z)); avi.l-z]
#    u = [fill(Inf, length(J[2])); zeros(length(J[3])+length(J[4])); zeros(length(z)); avi.u-z]
#
#    meaningful = find_non_trivial(A,l,u)
#    P2 = Poly([A[meaningful,:] spzeros(sum(meaningful), length(w))], l[meaningful], u[meaningful])
#
#    P2
#end
#
#function find_non_trivial(A,l,u)
#    [!(isinf(l[i]) && isinf(u[i])) || all(iszero.(A[i,:])) for i in 1:length(l)]
#end
#
#"""
#K[1] : Mz+Nw+o ≥ 0, z = l
#K[2] : Mz+Nw+o = 0, l ≤ z ≤ u
#K[3] : Mz+Nw+o ≤ 0, z = u
#K[4] : Mz+Nw+o free, l = z = u
#
## TODO probably inefficient initialization with Poly
#"""
#function local_piece(avi, n, m, K)
#    A = [avi.M avi.N;
#         I(n) spzeros(n,m)]
#
#    bounds = mapreduce(vcat, 1:n) do i
#        if i ∈ K[1]
#            [-avi.o[i] Inf avi.l[i] avi.l[i]]
#        elseif i ∈ K[2]
#            [-avi.o[i] -avi.o[i] avi.l[i] avi.u[i]] 
#        elseif i ∈ K[3]
#            [-Inf -avi.o[i] avi.u[i] avi.u[i]]
#        else
#            [-Inf Inf avi.l[i] avi.u[i]]
#        end
#    end
#    l = [bounds[:,1]; bounds[:,3]]
#    u = [bounds[:,2]; bounds[:,4]]
#    meaningful = find_non_trivial(A,l,u)
#    Poly(A[meaningful,:], l[meaningful], u[meaningful])
#end
#
#"""
#Form dictionary of index sets:
#J[1] = {i : lᵢ = zᵢ     , (Mz+Nw+o)ᵢ > 0 }
#J[2] = {i : lᵢ = zᵢ     , (Mz+Nw+o)ᵢ = 0 }
#J[3] = {i : lᵢ < zᵢ < uᵢ, (Mz+Nw+o)ᵢ = 0 }
#J[4] = {i :      zᵢ = uᵢ, (Mz+Nw+o)ᵢ = 0 }
#J[5] = {i :      zᵢ = uᵢ, (Mz+Nw+o)ᵢ < 0 }
#J[6] = {i : lᵢ = zᵢ = uᵢ, (Mz+Nw+o)ᵢ = 0 }
#"""
#function comp_indices(avi, z, w; tol=1e-4)
#    J = Dict{Int, Vector{Int}}()
#    r = avi.M*z+avi.N*w+avi.o
#    equal_bounds = isapprox.(avi.l, avi.u; atol=tol)
#    riszero = isapprox.(r, 0; atol=tol)
#    J[1] = findall( isapprox.(z, avi.l; atol=tol) .&& r .> tol )
#    J[2] = findall( isapprox.(z, avi.l; atol=tol) .&& riszero .&& .!equal_bounds)
#    J[3] = findall( (avi.l.+tol .< z .< avi.u.-tol) .&& riszero )
#    J[4] = findall( isapprox.(z, avi.u; atol=tol) .&& riszero .&& .!equal_bounds)
#    J[5] = findall( isapprox.(z, avi.u; atol=tol) .&& r .< -tol )
#    J[6] = findall( equal_bounds .&& riszero )
#    return J
#end
