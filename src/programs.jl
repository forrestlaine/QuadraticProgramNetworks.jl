struct Linear <: Function
    a::Vector{Float64} 
end

function Base.isequal(L1::Linear, L2::Linear)
    isequal(L1.a, L2.a)
end
function Base.hash(L::Linear, h::UInt)
    hash(("Linear",L.a), h)
end

function (f::Linear)(x::Vector{Float64})
    f.a'*x
end

struct Quadratic <: Function
    Q::SparseMatrixCSC{Float64, Int32}
    q::Vector{Float64}
end

function (f::Quadratic)(x::Vector{Float64})
    0.5*x'*(f.Q*x)+ x'*f.q
end

function Base.sum(fs::Union{Vector{Quadratic}, NTuple{N,Quadratic}}) where N
    Quadratic(sum(f.Q for f in fs), sum(f.q for f in fs)) 
end

struct QP
    f::Quadratic
    constraint_indices::Vector{Int}
    var_indices::Vector{Int}
end

"""
poly = S
group_mapping = {1 => 1, 2 => 1, 3 => 2}

implies that players 1 & 2 share a multiplier for constraint S,
and player 3 has their own multiplier. No other players are involved.
"""
struct Constraint
    poly::Poly
    group_mapping::Dict{Int, Int}
end

struct QEP
    qps::Dict{Int, QP}
    constraints::Dict{Int, Constraint}
end

@enum SharedVariableMode begin
    MIN_NORM = 1
    SHARED_DUAL = 2
end

Base.@kwdef mutable struct QPNetOptions
    shared_variable_mode::SharedVariableMode=SHARED_DUAL
    max_iters::Int=150
    tol::Float64=1e-4
    high_dimension::Bool=false
    high_dimension_max_iters::Int=10
    make_requests::Bool=false
    try_hull::Bool=true
    debug::Bool=true
    gen_solution_map::Bool=false
end

struct QPNet
    qps::Dict{Int, QP}
    constraints::Dict{Int, Constraint}
    network_edges::Dict{Int, Set{Int}}
    reachable_nodes::Dict{Int, Set{Int}}
    network_depth_map::Dict{Int, Set{Int}}
    options::QPNetOptions
    variables::Vector{Num}
    var_indices::Dict{Num, Int}
end
    
function QPNet(sym_vars::Vararg{Union{Num,Array{Num}}})
    all_vars = Num[] 
    var_inds = Dict{Num, Int}()
    ind = 0
    for sym_var in sym_vars
        for (i,var) in enumerate(sym_var)
            push!(all_vars, var)
            var_inds[var] = ind+i
        end
        ind += length(sym_var)
    end
    qps = Dict{Int, QP}()
    constraints = Dict{Int, Constraint}()
    network_edges =  Dict{Int, Set{Int}}()
    reachable_nodes =  Dict{Int, Set{Int}}()
    network_depth_map =  Dict{Int, Set{Int}}()
    options = QPNetOptions()
    QPNet(qps, constraints, network_edges, reachable_nodes, network_depth_map, options, all_vars, var_inds)    
end

function flatten(qpn::QPNet)
    QPNet(qpn.qps,
          qpn.constraints,
          Dict(1=>Set(keys(qpn.qps))),
          qpn.options,
          qpn.variables,
          qpn.var_indices)
end

function get_flat_initialization(qpn::QPNet; x0 = zeros(length(qpn.variables)))
    qpn_flat = flatten(qpn)
    qpn_flat.options.gen_solution_map = false
    ret = solve(qpn_flat, x0)
    ret.x_opt
end

"""
Convenience wrapper for setting up example problems. 
Usually not a good idea to abuse value types like this but
okay since these are not time-critical calls.
"""
function setup(sym::Symbol; args...)
    setup(Val(sym); args...)
end

variables(name, dims...) = Symbolics.variables(name, dims...)
variable(name) = Symbolics.variable(name)


function add_constraint!(qp_net, cons, lb, ub; tol=1e-8)
    @assert length(cons) == length(lb) == length(ub)
    A = Symbolics.sparsejacobian(cons, qp_net.variables)
    rows,cols,vals = findnz(A)
    n, m = size(A)
    try
        A = sparse(rows, cols, [Float64(v.val) for v in vals], n, m)
    catch err
        error("Detected non-linear constraint!")
    end
    droptol!(A, tol)

    vals = map(cons) do con
        vars = Symbolics.get_variables(con)
        Symbolics.substitute(con, Dict(var=>0.0 for var in vars)).val
    end
    
    poly = Poly(A,lb-vals,ub-vals)
    mapping = Dict{Int, Int}()
    constraint = Constraint(poly, mapping)
    id = maximum(keys(qp_net.constraints), init=0) + 1
    qp_net.constraints[id] = constraint
    return id
end

function add_qp!(qp_net, cost, con_inds, private_vars...; tol=1e-8)
    grad = Symbolics.gradient(cost, qp_net.variables)
    Q = Symbolics.sparsejacobian(grad, qp_net.variables)
    rows,cols,vals = findnz(Q)
    n,m = size(Q)
    try
        Q = sparse(rows, cols, [Float64(v.val) for v in vals], n, m)
    catch err
        error("Detected non-quadratic cost!")
    end
    q = map(x->Float64(x.val), Symbolics.substitute(grad, Dict(v => 0.0 for v in qp_net.variables)))
    droptol!(Q, tol)
    f = Quadratic(Q, q)
    var_inds = mapreduce(vcat, private_vars; init = Int[]) do var
        inds = Vector{Int}()
        foreach(vi->push!(inds, qp_net.var_indices[vi]), var)
        inds
    end
    player_id = maximum(keys(qp_net.qps), init=0) + 1
    qp_net.qps[player_id] = QP(f, con_inds, var_inds)
    #if level in keys(qp_net.network)
    #    push!(qp_net.network[level], player_id)
    #else
    #    qp_net.network[level] = Set(player_id)
    #end
    return player_id
end

"""
Deletes edges that are redundant in the adjacency matrix. Returns 
the minimal adjacency matrix as well as the reachability matrix,
where the R[i,j] = true if node j is reachable from node i.

If cyclic graphs are detected, returns an error, since cycles result in degenerate 
equilibrium points.

This is not the fastest implementation of this algorithm, but this is
for a non-critical path (setup).
"""
function create_minimal_adj_matrix(N, edge_list)
    A = zeros(Bool, N, N)
    for (i,j) in edge_list
        if i == j
            error("Cannot have self edges. (In this case, node $i -> $i).")
        end
        A[i,j] = true
    end

    R = zeros(Bool, N, N) # transition matrix
    Aⁿ = copy(A)
    edge_deleted = false
    for n = 2:N
        R .|= Aⁿ
        Aⁿ = (Aⁿ * A) .> 0
        for i = 1:N
            if Aⁿ[i,i]
                error("Cycle detected. (In this case, cycle leading from node $i -> $i after $n transitions.")
            end
            for j = 1:N
                if A[i,j] && Aⁿ[i,j]
                    A[i,j] = false
                    @info "Deleting $i -> $j (found alternate path using $n transitions)"
                end
            end
        end
    end
    A, R
end

"""
Uses the reachability matrix R to compute a depth map.
The nodes at depth d only have dependences on nodes of depth d+1 or greater. 
Any node in the graph can be reached from at least one of the nodes at depth 1.
"""
function create_depth_map(R)
    depth_map = Dict{Int, Set{Int}}()
    N = size(R,1)
    deleted_nodes = Set()
    d = 0
    Rd = copy(R)
    while length(deleted_nodes) < N
        nodes_at_depth = setdiff(Set([i for i in 1:N if !any(Rd[:,i])]), deleted_nodes)
        if length(nodes_at_depth) > 0
            d += 1
            depth_map[d] = nodes_at_depth
            union!(deleted_nodes, nodes_at_depth)
            remaining_nodes = setdiff(1:N, deleted_nodes)
            Rd = R[remaining_nodes,:]
        else
            error("Something appears wrong with the graph structure. There should always be nodes found for every increasing depth")
        end
    end
    @assert sum(sum(R[i,:] for i in depth_map[1]) .> 0) == N - length(depth_map[1])
    depth_map
end

"""
Add edges to the qp_net
"""
function add_edges!(qp_net, edge_list)
    N = length(qp_net.qps)
    A, R = create_minimal_adj_matrix(N, edge_list)
    depth_map = create_depth_map(R)
    for (d, nodes) in depth_map
        qp_net.network_depth_map[d] = nodes
    end
    for i in 1:N
        qp_net.network_edges[i] = Set(collect(1:N)[A[i,:]])
        qp_net.reachable_nodes[i] = Set(collect(1:N)[R[i,:]])
    end
end

"""
If group_map[constraint_id][player_a] == group_map[constraint_id][player_b], then 
player_a and player_b share a multiplier for constraint_id.

If no group_map is provided for a particular constraint, it is assumed that each player has a unique multiplier associated with the constriant.
"""
function assign_constraint_groups!(qp_net; group_map=Dict{Int, Dict{Int, Int}}())
    # TODO make dirty flag, set until this is called
    for (con_id, constraint) in qp_net.constraints
        for (player_id, qp) in qp_net.qps
            if con_id in qp.constraint_indices
                if con_id ∈ keys(group_map)
                    if player_id ∉ keys(group_map[con_id])
                        error("A group map is provided for constraint $con_id, but it does not contain mappings for all players who respect this constraint (namely player $player_id).")
                    end
                    group_id = group_map[con_id][player_id]
                else
                    group_id = player_id
                end
                constraint.group_mapping[player_id] = group_id
            end
        end
    end
end

function set_options!(qp_net; kwargs...)
    for kwarg in kwargs
        try
            setfield!(qp_net.options, kwarg[1], kwarg[2])
        catch e
            @warn("Invalid option name $(kwarg[1]) with value $(kwarg[2]), skipping")
        end
    end
end

function display_solution(qpn::QPNet, x)
    for var in qpn.variables
        idx = qpn.var_indices[var]
        val = x[idx]
        @info "($idx) $var => $val"
    end
end

function num_levels(qpn::QPNet)
    length(qpn.network_depth_map)
end

function gather(qpn::QPNet, level)
    qps = Dict(i=>qpn.qps[i] for i in qpn.network_depth_map[level])
    constraints = Dict{Int, Constraint}(id=>qpn.constraints[id] for qp in values(qps) for id in qp.constraint_indices)
    QEP(qps, constraints)
end

function fair_obj(qep::QEP)
    sum([qp.f for qp in values(qep.qps)])
end

function fair_obj(qpn::QPNet, level)
    sum([qpn.qps[i].f for i in qpn.network_depth_map[level]])
end

function level_indices(qpn::QPNet, level)
    reduce(vcat, (qpn.qps[i].var_indices for i in qpn.network_depth_map[level]))
end

function sub_indices(qpn::QPNet, level)
    L = length(qpn.network_depth_map)
    reduce(vcat, (qpn.qps[i].var_indices for l in level+1:L for i in qpn.network_depth_map[l]))
end

function subeq_indices(qpn::QPNet, level)
    L = length(qpn.network)
    reduce(vcat, (qpn.qps[i].var_indices for l in level:L for i in qpn.network[l]))
end

function param_indices(qpn::QPNet, level)
    collect(setdiff(1:embedded_dim(first(qpn.constraints).second.poly), Set(subeq_indices(qpn, level))))
end
