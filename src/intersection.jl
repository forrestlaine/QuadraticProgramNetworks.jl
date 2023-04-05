mutable struct IntersectionNode
    pus::Vector{PolyUnion}
    contributing_poly::Poly
    resulting_poly::Union{Poly, Missing}
    children::Vector{IntersectionNode}
    width::Int
    depth::Int
    empty::Bool
    level_dim::Int
    IntersectionNode(pus::Vector{PolyUnion}, depth, level_dim, state) = begin
        width = depth > length(pus) ? 0 : length(pus[depth])
        new(pus, pus[depth-1][state[depth-1]], missing, IntersectionNode[], width, depth, false, level_dim)
    end
end

mutable struct IntersectionRoot
    pus::Vector{PolyUnion}
    pms::Vector{Progress}
    show_progress::Bool
    children::Vector{IntersectionNode}
    width::Int
    leaf_count::Int
    depth_widths::Vector{Int}  # width i is total number of leaf nodes at or below any node at corresponding level 
    len::Int
    red_lengths::Vector{Int}
    level_dim::Int
    guide::Function
end

function Base.IteratorSize(root::IntersectionRoot)
    Base.SizeUnknown()
end
function Base.eltype(root::IntersectionRoot)
    Poly
end

function potential_length(root::IntersectionRoot)
    root.len
end

function potential_length(pu::PolyUnion)
    length(pu)
end

function depth(root::IntersectionRoot)
    length(root.pus)
end
depth(pu::PolyUnion) = 1

function set_guide!(root::IntersectionRoot, guide)
    #todo -- is there an intelligent way to guide intersection iteration?
    return
end

function IntersectionRoot(pus::Vector{PolyUnion}, red_lengths, level_dim; show_progress=false)
    @assert(length(pus) > 1)
    N = length(pus)
    pu_lengths = length.(pus)
    poly_inds = [ 1:L for L in pu_lengths ]
    depth_widths = [prod(pu_lengths[i+1:end]; init=1) for i in 1:N]
    children = IntersectionNode[]
    pms = [Progress(length(pu); desc="Depth "*string(i)*": ", dt=0.1, offset=i) for (i,pu) ∈ enumerate(pus)]
    IntersectionRoot(pus, pms, show_progress, children, pu_lengths[1], 0, depth_widths, prod(pu_lengths), red_lengths, level_dim, x->Inf)
end

function get_next!(node::IntersectionNode, parent_poly, state)
    if ismissing(node.resulting_poly)
        if isnothing(parent_poly)
            node.resulting_poly = deepcopy(node.contributing_poly)
        else
            node.resulting_poly = poly_intersect(node.contributing_poly, parent_poly)
        end
        if isempty(node.resulting_poly)
            node.empty = true
        end
    end
    if node.empty
        return nothing
    end
    if node.depth > length(state)
        if intrinsic_dim(node.resulting_poly) ≥ node.level_dim
            return node.resulting_poly 
        else
            node.empty = true
            return nothing
        end
    else
        while state[node.depth] ≤ node.width
            if state[node.depth] > length(node.children)
                push!(node.children, IntersectionNode(node.pus, node.depth+1, node.level_dim, state))
            end
            poly = get_next!(node.children[state[node.depth]], node.resulting_poly, state)
            if isnothing(poly)
                state[node.depth] += 1
                state[node.depth+1:end] .= 1
                continue
            else
                return poly
            end
        end
        return nothing # only if node is exhausted
    end
end

function Base.iterate(root::IntersectionRoot)
    if root.show_progress
        for i ∈ 1:length(root.pms)
            ProgressMeter.update!(root.pms[i], 1)
        end
    end
    iterate(root, ones(Int, length(root.pus)))
end

function Base.iterate(root::IntersectionRoot, state)
    while state[1] ≤ root.width
        if state[1] > length(root.children)
            push!(root.children, IntersectionNode(root.pus, 2, root.level_dim, state))
        end
        poly = get_next!(root.children[state[1]], nothing, state)
        full_lengths = length.(root.pus)
        redzone = all(ind > full_len-red_len for (ind, red_len, full_len) in zip(state, root.red_lengths, full_lengths))
        root.leaf_count = sum((ind-1)*d for (ind,d) in zip(state, root.depth_widths))
        increment!(state, length.(root.pus))
        if root.show_progress
            for i ∈ 1:length(root.pms)
                ProgressMeter.update!(root.pms[i], state[i])
            end
        end
        if isnothing(poly) || redzone # this should happen only when child is exhausted
            continue
        end
        return (poly, state)
    end

    return nothing
end

function increment!(state, lens)
    N = length(lens)
    for n = N:-1:1
        offset = n==1 ? 1 : 0
        if state[n] < lens[n] + offset
            state[n] += 1
            return
        else
            state[n] = 1
        end
    end
end
