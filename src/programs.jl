struct Quadratic <: Function
    Q::SparseMatrixCSC{Float64, Int32}
    q::Vector{Float64}
end

function (f::Quadratic)(x::Vector{Float64})
    0.5*x'*(f.Q*x + f.q)
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
    debug::Bool=true
    gen_solution_map::Bool=false
end

struct QPNet
    qps::Dict{Int, QP}
    constraints::Dict{Int, Constraint}
    network::Dict{Int, Set{Int}}
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
    network =  Dict{Int, Set{Int}}()
    options = QPNetOptions()
    QPNet(qps, constraints, network, options, all_vars, var_inds)    
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
    ret = solve(qpn_flat, x0)
    ret.x_opt
end

function variables(name, dims...)
    Symbolics.variables(name, dims...)
end


function add_constraint!(qp_net, cons, lb, ub; tol=1e-8)
    @assert length(cons) == length(lb) == length(ub)
    A = Symbolics.sparsejacobian(cons, qp_net.variables)
    rows,cols,vals = findnz(A)
    n, m = size(A)
    A = sparse(rows, cols, [Float64(v.val) for v in vals], n, m)
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

function add_qp!(qp_net, level, cost, con_inds, private_vars...; tol=1e-8)
    grad = Symbolics.gradient(cost, qp_net.variables)
    Q = Symbolics.sparsejacobian(grad, qp_net.variables)
    rows,cols,vals = findnz(Q)
    n,m = size(Q)
    Q = sparse(rows, cols, [Float64(v.val) for v in vals], n, m)
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
    if level in keys(qp_net.network)
        push!(qp_net.network[level], player_id)
    else
        qp_net.network[level] = Set(player_id)
    end
    return player_id
end

function assign_constraint_groups!(qp_net; group_map=Dict{Int, Dict{Int, Int}}())
    # TODO make dirty flag, set until this is called
    for (con_id, constraint) in qp_net.constraints
        for (player_id, qp) in qp_net.qps
            if con_id in qp.constraint_indices
                group_id = 
                    (con_id ∈ keys(group_map)) && (player_id ∈ keys(group_map[con_id])) ? 
                    group_map[con_id][player_id] : 
                    player_id
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



function num_levels(qpn::QPNet)
    length(qpn.network)
end

function gather(qpn::QPNet, level)
    qps = Dict(i=>qpn.qps[i] for i in qpn.network[level])
    constraints = Dict{Int, Constraint}(id=>qpn.constraints[id] for qp in values(qps) for id in qp.constraint_indices)
    QEP(qps, constraints)
end

function fair_obj(qep::QEP)
    sum([qp.f for qp in values(qep.qps)])
end

function fair_obj(qpn::QPNet, level)
    sum([qpn.qps[i].f for i in qpn.network[level]])
end

function sub_indices(qpn::QPNet, level)
    L = length(qpn.network)
    reduce(vcat, (qpn.qps[i].var_indices for l in level+1:L for i in qpn.network[l]))
end

function subeq_indices(qpn::QPNet, level)
    L = length(qpn.network)
    reduce(vcat, (qpn.qps[i].var_indices for l in level:L for i in qpn.network[l]))
end

function param_indices(qpn::QPNet, level)
    collect(setdiff(1:embedded_dim(first(qpn.constraints).second.poly), Set(subeq_indices(qpn, level))))
end
