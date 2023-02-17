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
Represents a generalized affine variational inequality
(split into two conditions for computational reasons)

(Mz + Nw + o) ⟂ (l₁ ≤  z₁   ≤ u₁)
(     z₂    ) ⟂ (l₂ ≤ Az+Bw ≤ u₂)

Possible todo: add support for following conditions.
(M₃z + N₃w + o₃) ⟂ (l₃ ≤ A₃z+B₃w ≤ u₃)
"""
struct GAVI
    M::SparseMatrixCSC{Float64, Int32}
    N::SparseMatrixCSC{Float64, Int32}
    o::Vector{Float64}
    l1::Vector{Float64}
    u1::Vector{Float64}
    A::SparseMatrixCSC{Float64, Int32}
    B::SparseMatrixCSC{Float64, Int32}
    l2::Vector{Float64}
    u2::Vector{Float64}
end

"""
Given M, L, o, l, u, w,
Find z, u, v, s.t.:
    u - v = M z + N w + o
    u ≥ 0 ⟂ z - l ≥ 0
    v ≥ 0 ⟂ u - z ≥ 0
Currently uses PATHSolver
"""
function solve_avi(avi::AVI, z0, w)
    PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
    (path_status, z, info) =  PATHSolver.solve_mcp(avi.M, avi.N*w+avi.o,avi.l, avi.u, z0, silent=true, convergence_tolerance=1e-8)
    (; sol_bad, degree, r) = check_avi_solution(avi, z, w)
    @infiltrate sol_bad
    status = (path_status == PATHSolver.MCP_Solved || path_status == PATHSolver.MCP_Solved) ? SUCCESS : FAILURE
    status == FAILURE && @infiltrate
    return (; z, status)
end

function solve_avi(gavi::GAVI, z0, w)
    avi = convert(gavi)
    d2 = length(gavi.l2)
    z0 = [z0; zeros(d2)]
    (; z, status) = solve_avi(avi, z0, w)
    z = z[1:end-d2]
    (; z, status)
end

function convert(gavi::GAVI)
    d1 = length(gavi.l1)
    d2 = length(gavi.l2)
    M = [gavi.M spzeros(d1,d2);
         gavi.A -sparse(I,d2,d2);
         spzeros(d2,d1) sparse(I, d2,d2) spzeros(d2,d2)]
    N = [gavi.N; gavi.B; spzeros(d2, size(gavi.N,2))]
    o = [gavi.o; zeros(d2); zeros(d2)]
    l = [gavi.l1; fill(-Inf, d2); gavi.l2]
    u = [gavi.u1; fill(Inf, d2); gavi.u2]
    AVI(M,N,o,l,u)
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
function solve_qep(qep_base, x, S=nothing, shared_decision_inds=Vector{Int}(); 
                   debug=false,
                   high_dimension=false,
                   shared_variable_mode=SHARED_DUAL,
                   rng=MersenneTwister(1))

    x_dim = length(x)
    N_players = length(qep_base.qps)

    if isnothing(S)
        qep = qep_base
        aux_dim = 0
        @assert length(shared_decision_inds) == 0
    else
        qep = deepcopy(qep_base)
        foreach(qp->push!(qp.constraint_indices, -1), values(qep.qps))
        if shared_variable_mode==SHARED_DUAL
            qep.constraints[-1] = Constraint(S, Dict(i=>1 for i in keys(qep.qps))) 
        else
            qep.constraints[-1] = Constraint(S, Dict(i=>i for i in keys(qep.qps))) 
        end
        aux_dim = embedded_dim(S) - x_dim
        @assert length(shared_decision_inds) > 0
    end

    private_decision_inds = reduce(vcat, (qp.var_indices for qp in values(qep.qps))) |> sort
    decision_inds = [private_decision_inds; shared_decision_inds]
    
    N_shared_vars = length(shared_decision_inds) + aux_dim
    N_private_vars = sum(length.(private_decision_inds))
   
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
        group_labels = keys(group_to_player_map)
        num_groups = length(group_labels)

        A1 = repeat(A[:, decision_inds], num_groups)
        Ax = (local_aux_dim > 0) ? repeat(A[:, x_dim+1:x_dim+local_aux_dim], num_groups) : spzeros(num_groups*length(l), aux_dim)
        B1 = repeat(A[:, param_inds], num_groups)
        l1 = repeat(l, num_groups)
        u1 = repeat(u, num_groups)
        A2 = mapreduce(vcat, group_labels) do gid
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

    Ma = [M11 M12 M13 M14;
         M21 M22 M23 M24;
         M31 M32 M33 M34;
         M41 M42 M43 M44]

    M = [M11 M12 M14;
         M21 M22 M24]

    N1 = Rs
    N2 = spzeros(N_shared_vars, length(param_inds))
    N3 = spzeros(standard_dual_dim, length(param_inds))
    N4 = Bs

    Na = [N1; N2; N3; N4]

    N = [N1; N2]

    oa = [qs; zeros(N_shared_vars+2*standard_dual_dim)]

    o = [qs; zeros(N_shared_vars)]

    l = [fill(-Inf, length(qs)+N_shared_vars); ls; fill(-Inf, standard_dual_dim)]
    u = [fill(Inf, length(qs)+N_shared_vars); us; fill(Inf, standard_dual_dim)]

    l1 = fill(-Inf, length(qs) + N_shared_vars)
    u1 = fill(Inf, length(qs) + N_shared_vars)

    A = [M41 M42 M44]
    B = N4
    l2 = ls
    u2 = us


    w = x[param_inds]
    z0a = [x[decision_inds]; zeros(aux_dim+N_players*N_shared_vars); M41*[x[decision_inds]; zeros(aux_dim)]+N4*w; zeros(standard_dual_dim)]

    z0 = [x[decision_inds]; zeros(aux_dim+N_players*N_shared_vars); zeros(standard_dual_dim)]
    gavi = GAVI(M,N,o,l1,u1,A,B,l2,u2)
    avi = AVI(Ma, Na, oa, l, u)
    (; z, status) = solve_avi(gavi, z0, w)
    status != SUCCESS && @infiltrate
    status != SUCCESS && error("AVI solve error!")
    ψ_inds = collect(N_private_vars+N_shared_vars+1:N_private_vars+N_shared_vars*(N_players+1))

    if shared_variable_mode == MIN_NORM
        if high_dimension 
            (; piece, x_opt, reduced_inds) = get_single_avi_solution(gavi,z,w,decision_inds,param_inds,rng; debug)
            if length(ψ_inds) > 0
                old_piece = piece
                old_x_opt = x_opt
                ## IMPORTANT: not passing in ψ_inds, since by default ψ_inds will
                ## correspond to the first indices immediately after all x_dims in piece
                num_remaining_ψ = length(ψ_inds) - length(intersect(ψ_inds, reduced_inds)) # ψ_inds and reduced_inds are both pre-permutation
                f_min_norm = min_norm_objective(embedded_dim(piece), length(x_opt), num_remaining_ψ)
                (; piece, x_opt, z_revised) = revise_avi_solution(f_min_norm, piece, x_opt, decision_inds, param_inds, rng)
                println("ψ_val: ", norm(z_revised[length(x_opt)+1:length(x_opt)+num_remaining_ψ]))
                new_dim = embedded_dim(piece)
                old_dim = embedded_dim(old_piece)
                println("Piece increased in size to ", new_dim, " from ", old_dim)
            end
            (; x_opt, Sol=[piece,], f_up=nothing)
        else
            @error "not implemented yet" 
        end
    else
        if high_dimension
            (; piece, x_opt) = get_single_avi_solution(avi,z,w,decision_inds,param_inds,rng; debug)
            (; x_opt, Sol=[piece,], f_up=nothing)
        else 
            x_opt = copy(x)
            x_opt[decision_inds] = z[1:length(decision_inds)]
            Sol = LocalAVISolutions(avi, z, w, decision_inds, param_inds)
            (; x_opt, Sol, f_up=nothing)
        end
    end
end
"""
f(z) = 0.5*∑(zᵢ²; i∈1...n if m<i≤m+l)
"""
function min_norm_objective(n, m, l)
    Q = spzeros(n,n)
    Q[m+1:m+l, m+1:m+l] = sparse(I, l, l)
    Quadratic(Q, zeros(n))
end

function revise_avi_solution(f, piece, x, decision_inds, param_inds, rng)
    # TODO refactor this to use solve_qep (need to call this function from
    # algorithm.jl)
    # TODO update to use GAVI formulation (keep dims down)
    # TODO Need to make sure that x_param values don't change.
    (A, ll, uu, _, _) = vectorize(piece)
    (m,n) = size(A)

    non_param_inds = setdiff(Set(collect(1:n+m)), Set(param_inds)) |> collect |> sort
    npinds_1 = non_param_inds ∩ (1:n)
    npinds_2 = non_param_inds ∩ (n+1:n+m)
    J = [f.Q -A';
         A spzeros(m,m)]
    l = [fill(-Inf, n); ll]
    u = [fill(Inf, n); uu]
    M = J[npinds_1, non_param_inds]
    N = J[npinds_1, param_inds]
    A = J[npinds_2, non_param_inds]
    B = J[npinds_2, param_inds]

    o = f.q[npinds_1]
    l1 = l[npinds_1]
    u1 = u[npinds_1]
    l2 = l[npinds_2]
    u2 = u[npinds_2]
    
    gavi = GAVI(M, N, o, l1, u1, A, B, l2, u2)
    z0 = [x; zeros(length(f.q)-length(x)); zeros(m)][non_param_inds]

    (; z, status) = solve_avi(gavi, z0, x[param_inds])
    status != SUCCESS && @infiltrate
    status != SUCCESS && error("AVI solve error!")
    (; piece, x_opt, reduced_inds) = get_single_avi_solution(gavi, z, x[param_inds], sort(decision_inds), param_inds, rng)
    (; piece, x_opt, z_revised=z)
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
