mutable struct IntersectionRoot
    pus::Vector{PolyUnion2}
    children::Vector{IntersectionNode}
    width::Int
    leaf_count::Int
    depth_widths::Vector{Int}  # width i is total number of leaf nodes at or below any node at corresponding level 
    len::Int
    red_lengths::Vector{Int}
end

function Base.length(root::IntersectionRoot)
    root.len
end

mutable struct IntersectionNode
    pus::Vector{PolyUnion2}
    contributing_poly::Poly2
    resulting_poly::Union{Poly2, Missing}
    children::Vector{IntersectionNode}
    width::Int
    depth::Int
    used::Bool
    IntersectionNode(pus::Vector{PolyUnion2}, depth, state) = begin
        new(pus, pus[depth-1][state[depth-1]], missing, IntersectionNode[], length(pus[depth]), depth, false)
    end
end

function IntersectionRoot(pus::Vector{PolyUnion2}, red_lengths)
    N = length(pus)
    pu_lengths = length.(pus)
    poly_inds = [ 1:L for L in pu_lengths ]
    depth_widths = [prod(pu_lengths[i+1:end]; init=1) for i in 1:N]
    children = IntersectionNode[]
    IntersectionRoot(pus, children, 1, 0, depth_widths, prod(pu_lengths), red_lengths)
end


function get_next!(node::IntersectionNode, parent_poly, state)
    if ismissing(node.resulting_poly)
        if isnothing(parent_poly)
            node.resulting_poly = node.contributing_poly
        else
            node.resulting_poly = poly_intersect(node.contributing_poly, parent_poly)
        end
        if isempty(node.resulting_poly)
            return nothing
        end
    end
    if node.used
        return nothing
    end
    if node.depth == length(state)
        node.used = true
        return node.resulting_poly 
    else
        while state[node.depth] ≤ node.width
            if state[node.depth] > length(node.children)
                push!(node.children, IntersectionNode(node.pus, node.depth+1, state))
            end
            poly = get_next!(node.children[state[node.depth]], node.resulting_poly, state)
            state[node.depth] += 1
            isnothing(poly) && continue
            return poly
        end
        return nothing
    end
end

function Base.iterate(root::IntersectionRoot)
    iterate(root, ones(Int, length(root.pus)))
end

function Base.iterate(root::IntersectionRoot, state)
    while state[1] ≤ root.width
        if state[1] ≥ length(root.children)
            push!(root.children, IntersectionNode(root.pus, 2, state)
        end
        poly = get_next!(root.children[state[1]], nothing, state)
        @infiltrate
        if isnothing(poly)
            state[1] += 1
            continue
        else
            if all(ind ≤ red_len for (ind, red_len) in zip(state, root.red_lengths))
                increment!(state, length.(root.pus))
                continue
            end
            root.leaf_count = sum((ind-1)*d for (ind,d) in zip(state, root.depth_widths))
            increment!(state, length.(root.pus))
            return (node, state)
        end
    end
    return nothing
end

function increment!(state, lens)
    N = length(lens)
    for n = N:-1:1
        if state[n] < lens[n]
            state[n] += 1
            return
        else
            state[n] = 1
        end
    end
    error("Shouldn't ever loop back around like this.")
end




##########################

mutable struct IntersectionRoot
    children
    child_index::Int
    leaf_count::Int
    depth_widths::Vector{Int}  # width i is total number of leaf nodes at or below any node at corresponding level 
    len::Int
    red_lengths::Vector{Int}
end

function Base.length(root::IntersectionRoot)
    root.len
end

mutable struct IntersectionNode
    parent::Union{IntersectionNode, IntersectionRoot}
    contributing_poly::Poly2
    resulting_poly::Union{Poly2, Missing}
    children::Vector{IntersectionNode}
    child_index::Int
    used::Bool
    IntersectionNode(poly, pus::Vector{PolyUnion2}) = begin
        node = new()
        node.contributing_poly = poly
        node.resulting_poly = missing
        pu_lengths = length.(pus)
        if length(pus) == 0
            children = IntersectionNode[]
        else
            children = [IntersectionNode(pus[1][i], pus[2:end]) for i in 1:pu_lengths[1]]
        end
        node.children = children
        node.child_index = 1
        node.used = false
        node
    end
end

function IntersectionRoot(pus::Vector{PolyUnion2}, red_lengths)
    N = length(pus)
    pu_lengths = length.(pus)
    poly_inds = [ 1:L for L in pu_lengths ]
    depth_widths = [prod(pu_lengths[i+1:end]; init=1) for i in 1:N]
    @infiltrate
    children = [IntersectionNode(pus[1][i], pus[2:end]) for i in 1:pu_lengths[1]]
    IntersectionRoot(children, 1, 0, depth_widths, prod(pu_lengths), red_lengths) |> assign_parent_to_children!
end

function assign_parent_to_children!(parent)
    for child in parent.children
        child.parent = parent
        assign_parent_to_children!(child)
    end
    parent
end

function get_next!(node::IntersectionNode)
    if ismissing(node.resulting_poly)
        if node.parent isa IntersectionRoot
            node.resulting_poly = node.contributing_poly
        else
            node.resulting_poly = poly_intersect(node.contributing_poly, node.parent.resulting_poly)
        end
        if isempty(node.resulting_poly)
            return nothing
        end
    end
    if node.used
        return nothing
    end
    if length(node.children) == 0
        node.used = true
        return (node.resulting_poly, [node.parent.child_index,])
    else
        while node.child_index ≤ length(node.children)
            child_result = get_next!(node.children[node.child_index])
            node.child_index += 1
            if isnothing(child_result)
                continue
            else
                poly = child_result[1]
                inds = child_result[2]
                push!(inds, node.parent.child_index)
                return (poly, inds)
            end
        end
        return nothing
    end
end

function Base.iterate(root::IntersectionRoot)
    while root.child_index ≤ length(root.children)
        result = get_next!(root.children[root.child_index])
        if isnothing(result)
            root.child_index += 1
            continue
        else
            (node, indices) = result
            if all(ind ≤ red_len for (ind, red_len) in zip(indices, root.red_lengths))
                return iterate(root)
            else
                root.leaf_count = sum((ind-1)*d for (ind,d) in zip(indices, root.depth_widths))
                return (node, nothing)
            end
        end
    end
    return nothing
end
function Base.iterate(root::IntersectionRoot, state)
    iterate(root)
end
