@enum StatusCode begin
    SUCCESS = 1
    RAY_TERM = 2
    MAX_ITERS = 3
    FAILURE = 4
end

ew_vcat = (args...) -> vcat.(args...)

struct AVI
    M::SparseMatrixCSC{Float64, Int32}
    N::SparseMatrixCSC{Float64, Int32}
    o::Vector{Float64}
    l::Vector{Float64}
    u::Vector{Float64}
end

"""
Given M, L, o, l, u, w,
Find z, u, v, s.t.:
    u - v = M z + N w + o
    u ≥ 0 ⟂ z - l ≥ 0
    v ≥ 0 ⟂ u - z ≥ 0
Currently uses PATHSolver
"""
function solve_avi(avi, z0, w)
    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    (path_status, z, info) =  PATHSolver.solve_mcp(avi.M, avi.N*w+avi.o,avi.l, avi.u, z0, silent=false, convergence_tolerance=1e-8)
    (; sol_bad, degree, r) = check_avi_solution(avi, z, w)
    @infiltrate sol_bad
    status = (path_status == PATHSolver.MCP_Solved || path_status == PATHSolver.MCP_Solved) ? SUCCESS : FAILURE
    status == FAILURE && @infiltrate
    return (; z, status)
end

function check_avi_solution(avi, z, w; tol=1e-6)
    r = avi.M*z + avi.N*w + avi.o
    r_pos = r .> tol
    r_neg = r .< -tol
    bad_count = sum(abs.(z[r_pos]-avi.l[r_pos]) .> tol) + sum(abs.(z[r_neg]-avi.u[r_neg]) .> tol)
    return (; sol_bad = bad_count > 0, degree = bad_count, r)
end


"""
Solve the Quadratic Equilibrium Problem.
"""
function solve_qep(qep_base, x, S=nothing, shared_decision_inds=Vector{Int}(); debug=false, high_dimension=false, shared_var_mode=SHARED_DUAL, rng=MersenneTwister(1)) 
    x_dim = length(x)
    N_players = length(qep_base.qps)

    if isnothing(S)
        qep = qep_base
        aux_dim = 0
        @assert length(shared_decision_inds) == 0
    else
        qep = deepcopy(qep_base)
        foreach(qp->push!(qp.constraint_indices, -1), values(qep.qps))
        qep.constraints[-1] = Constraint(S, Dict(i=>1 for i in keys(qep.qps))) 
        aux_dim = embedded_dim(S) - x_dim
        @assert length(shared_decision_inds) > 0
    end


    private_decision_inds = reduce(vcat, (qp.var_indices for qp in values(qep.qps))) |> sort
    decision_inds = [private_decision_inds; shared_decision_inds]
    
    N_shared_vars = length(shared_decision_inds) + aux_dim
    N_private_vars = sum(length.(private_decision_inds))
   
    # total_dual_dim = sum(length(S) for S in values(qep.sets); init=0)
    # TODO shared variable length needs to include auxiliary variables?
    standard_dual_dim = sum(length(Set(values(C.group_mapping)))*length(C.poly) for C in values(qep.constraints))

    param_inds = setdiff(Set(1:x_dim), Set(decision_inds)) |> collect |> sort

    player_order = collect(keys(qep.qps)) |> sort
    constraint_order = collect(keys(qep.constraints)) |> sort
    
    Qs, Mψ, Rs, qs = map(player_order) do i
        qp = qep.qps[i]
        inds = [qp.var_indices; shared_decision_inds]
        Q = [[qp.f.Q[inds, decision_inds]; spzeros(aux_dim, length(decision_inds))] spzeros(length(inds)+aux_dim, aux_dim)]
        M = mapreduce(hcat, player_order) do j
            j == i ? 
            [spzeros(length(qp.var_indices), N_shared_vars); -sparse(I, N_shared_vars, N_shared_vars)] :
            spzeros(length(qp.var_indices)+N_shared_vars, N_shared_vars)
        end
        R = [qp.f.Q[inds, param_inds]; spzeros(aux_dim, length(param_inds))]
        q = [qp.f.q[inds]; zeros(aux_dim)]
        Q, M, R, q
    end |> (x->vcat.(x...))
    

    As, A2s, Bs, ls, us = map(constraint_order) do id
        # TODO make sure that As / A2s accounts for auxiliary variables
        # properly (should be part of shared vars I think)
        A,l,u,_,_ = vectorize(qep.constraints[id].poly)
        local_aux_dim = size(A,2) - x_dim
        player_to_group_map = qep.constraints[id].group_mapping
        group_to_player_map = Dict{Int, Vector{Int}}()
        for (player, group) in player_to_group_map
            if group ∈ keys(group_to_player_map)
                push!(group_to_player_map[group], player)
            else
                group_to_player_map[group] = [player,]
            end
        end
        num_groups = length(group_to_player_map)

        A1 = repeat(A[:, decision_inds], num_groups)
        Ax = (local_aux_dim > 0) ? repeat(A[:, x_dim+1:x_dim+local_aux_dim], num_groups) : spzeros(num_groups*length(l), aux_dim)
        B1 = repeat(A[:, param_inds], num_groups)
        l1 = repeat(l, num_groups)
        u1 = repeat(u, num_groups)
        A2 = mapreduce(vcat, 1:num_groups) do gid
            mapreduce(hcat, player_order) do pid
                inds = [qep.qps[pid].var_indices; shared_decision_inds]
                if pid ∈ group_to_player_map[gid]
                    Ainds = A[:, inds]
                    Aaux = (local_aux_dim > 0) ? A[:, x_dim+1:x_dim+local_aux_dim] : spzeros(length(l), aux_dim)
                    [Ainds Aaux]
                else
                    spzeros(length(l), length(inds)+aux_dim)
                end
                #inds = [qep.qps[pid].var_indices; shared_decision_inds; x_dim+1:x_dim+local_aux_dim]
                #pid ∈ group_to_player_map[gid] ? A[:, inds] : spzeros(length(l), length(inds))
            end
        end
        #Ax = mapreduce(hcat, constraint_order) do i
        #    (i == id) ? A[:, x_dim+1:end] : spzeros(length(l), auxiliary_dims[i])
        #end |> (mat -> repeat(mat, num_groups, 1))

        #A2 = mapreduce(hcat, player_order) do i
        #    id ∈ keys(qep.qps[i].S) ? qep.qps[i].S[id] * (A[:, qep.qps[i].indices]) : spzeros(length(l), length(qep.qps[i].indices))
        #end
        #Ax = mapreduce(hcat, constraint_order) do i
        #    (i == id) ? A[:, x_dim+1:end] : spzeros(length(l), auxiliary_dims[i])
        #end
        [A1 Ax], A2, B1, l1, u1
    end |> (x->vcat.(x...))

    M11 = Qs
    M12 = Mψ
    M13 = spzeros(N_players*N_shared_vars+N_private_vars, standard_dual_dim)
    M14 = -A2s'
    M21 = spzeros(N_shared_vars, N_private_vars+N_shared_vars)
    M22 = repeat(sparse(I, N_shared_vars, N_shared_vars), 1, N_players)
    M23 = spzeros(N_shared_vars, standard_dual_dim)
    M24 = spzeros(N_shared_vars, standard_dual_dim)
    M31 = spzeros(standard_dual_dim, N_private_vars+N_shared_vars)
    M32 = spzeros(standard_dual_dim, N_players*N_shared_vars)
    M33 = spzeros(standard_dual_dim, standard_dual_dim)
    M34 = sparse(I, standard_dual_dim, standard_dual_dim)
    M41 = As
    M42 = spzeros(standard_dual_dim, N_shared_vars*N_players)
    M43 = -sparse(I, standard_dual_dim, standard_dual_dim)
    M44 = spzeros(standard_dual_dim, standard_dual_dim)

    M = [M11 M12 M13 M14;
         M21 M22 M23 M24;
         M31 M32 M33 M34;
         M41 M42 M43 M44]

    N1 = Rs
    N2 = spzeros(N_shared_vars, length(param_inds))
    N3 = spzeros(standard_dual_dim, length(param_inds))
    N4 = Bs

    N = [N1; N2; N3; N4]

    o = [qs; zeros(N_shared_vars+2*standard_dual_dim)]
    l = [fill(-Inf, length(qs)+N_shared_vars); ls; fill(-Inf, standard_dual_dim)]
    u = [fill(Inf, length(qs)+N_shared_vars); us; fill(Inf, standard_dual_dim)]

    w = x[param_inds]
    z0 = [x[decision_inds]; zeros(aux_dim+N_players*N_shared_vars); M41*[x[decision_inds]; zeros(aux_dim)]+N4*w; zeros(standard_dual_dim)]
    avi = AVI(M, N, o, l, u)

    (; z, status) = solve_avi(avi, z0, w)
    status != SUCCESS && @infiltrate
    status != SUCCESS && error("AVI Solve error!")
    if high_dimension
        (; piece, x_opt) = get_single_avi_solution(avi,z,w,decision_inds,param_inds,rng; debug)
        (; x_opt, Sol=[piece,])
    else 
        x_opt = copy(x)
        x_opt[decision_inds] = z[1:length(decision_inds)]
        Sol = LocalAVISolutions(avi, z, w, decision_inds, param_inds)
        (; x_opt, Sol)
    end
end

#function solve_qep(qep, x, S, sub_inds; debug=false, high_dimension=false, rng=MersenneTwister(1))
#    # TODO why not just preemptively add all subpiece sets to the dictionary,
#    # and only modify the qp constraint dependencies? This requires a ton of
#    # extra copying of huge matrices for large problems. Too lazy to fix now.
#    qep_augmented = deepcopy(qep)
#    foreach(qp->push!(qp.S, -1), values(qep_agumented.qps))
#    qep_agumented.sets[-1] = Constraint(S, Dict(i=>1 for i in keys(qep_augmented.qps)))
#
#    foreach(qp->qp.S[-1]=1.0, values(qep_augmented.qps))
#    qep_augmented.sets[-1] = S
#
#    qp_fair = QP(fair_obj(qep), Dict(id=>sum(qp.S[id] for qp in values(qep_augmented.qps) if id ∈ keys(qp.S); init=0.0) for id in keys(qep_augmented.sets)), sub_inds)
#    filter!(p->!iszero(p.second), qp_fair.S)
#    qep_augmented.qps[-1] = qp_fair
#    solve_qep(qep_augmented, x; debug, high_dimension)
#end
