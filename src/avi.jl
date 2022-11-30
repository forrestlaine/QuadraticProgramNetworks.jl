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
    (path_status, z, info) =  PATHSolver.solve_mcp(avi.M, avi.N*w+avi.o,avi.l, avi.u, z0, silent=true)
    status = (path_status == PATHSolver.MCP_Solved || path_status == PATHSolver.MCP_OK) ? SUCCESS : FAILURE
    status == FAILURE && @infiltrate
    return (; z, status)
end

function project_and_permute(S, var_inds, param_inds)
    d = embedded_dim(S)
    dv = length(var_inds)
    dp = length(param_inds)
    projection_inds = union(1:dv, d-dp+1:d)
    
    piece = project(S, projection_inds)
     
    #Atmp = similar(piece.A)
    #Atmp[:,var_inds] .= piece.A[:, 1:dv]
    #Atmp[:,param_inds] .= piece.A[:,dv+1:end]
    #piece.A .= Atmp
    #dropzeros!(piece.A)
    #meaningful = find_non_trivial(piece.A, piece.l, piece.u)
    #Poly(piece.A[meaningful,:], piece.l[meaningful], piece.u[meaningful])
    
    for slice in piece
        atmp = copy(slice.a)
        atmp[var_inds] = slice.a[1:dv]
        atmp[param_inds] = slice.a[dv+1:end]
        slice.a .= atmp
        dropzeros!(slice.a)
    end
    simplify(piece)
end

"""
Solve the Quadratic Equilibrium Problem.
"""
function solve_qep(qep, x; debug=false) 
    st = time() 
    decision_inds = reduce(vcat, (qp.indices for qp in values(qep.qps))) |> sort
    total_dual_dim = sum(length(S) for S in values(qep.sets); init=0)
    param_inds = setdiff(Set(1:length(first(values(qep.qps)).f.q)), Set(decision_inds)) |> collect |> sort

    player_order = collect(keys(qep.qps)) |> sort
    constraint_order = collect(keys(qep.sets)) |> sort

    Qs, Rs, qs = map(player_order) do i
        qp = qep.qps[i]
        qp.f.Q[qp.indices, decision_inds], qp.f.Q[qp.indices, param_inds], qp.f.q[qp.indices]
    end |> (x->vcat.(x...))
    
    As, A2s, Bs, ls, us = map(constraint_order) do id
        A,l,u,_,_ = vectorize(qep.sets[id])
        A2 = mapreduce(vcat, player_order) do i
            id ∈ keys(qep.qps[i].S) ? qep.qps[i].S[id] * (-A[:, qep.qps[i].indices])' : spzeros(length(qep.qps[i].indices), length(l))
        end
        A[:,decision_inds], A2', A[:, param_inds], l, u 
    end |> (x->vcat.(x...))

    M11 = Qs
    M12 = spzeros(length(decision_inds), total_dual_dim)
    M13 = A2s'
    #M21 = M12'
    M22 = spzeros(total_dual_dim, total_dual_dim)
    M23 = I(total_dual_dim)
    M31 = As
    #M32 = -M23'
    #M33 = M22

    M = [M11 M12 M13;
         M12' M22 M23;
         M31 -M23' M22]

    N1 = Rs
    N2 = spzeros(total_dual_dim, length(param_inds))
    N3 = Bs

    N = [N1; N2; N3]

    o = [qs; zeros(2*total_dual_dim)]
    l = [fill(-Inf, length(decision_inds)); ls; fill(-Inf, total_dual_dim)]
    u = [fill(Inf, length(decision_inds)); us; fill(Inf, total_dual_dim)]

    w = x[param_inds]
    z0 = [x[decision_inds]; zeros(total_dual_dim); M31*x[decision_inds]+N3*w]
    avi = AVI(M, N, o, l, u)
    en = time()
    debug && println("setting up avi took ", en-st)
    
    (; z, status) = solve_avi(avi, z0, w)
    x_opt = copy(x)
    x_opt[decision_inds] = z[1:length(decision_inds)]

    status != SUCCESS && @infiltrate
    status != SUCCESS && error("AVI Solve error!")
    #Sol, single_Sol = local_pieces(avi, z, w)
    st = time() 
    debug && println("solving avi took ", st-en)
    #xinds = [collect(1:length(decision_inds)); collect(length(z)+1:length(z)+length(w))]
    Sol = LocalAVISolutions(avi, z, w, decision_inds, param_inds)
    en = time()
    debug && println("Creating sol set took ", en-st)
    #local_avi_sols = LocalAVISolutions(avi, z, w, decision_inds, param_inds)
    #Sol = Iterators.filter(!isempty, (project_and_permute(S, decision_inds, param_inds) for S in local_avi_sols))
    
    #Sol = Iterators.filter(!isempty, (project_and_permute(S, all_inds, param_inds) for S in local_pieces(avi, z, w; debug, expansion)))
    #single_Sol = project_and_permute(single_Sol, all_inds, param_inds)
    
    x_opt = copy(x)
    x_opt[decision_inds] = z[1:length(decision_inds)]
    (; x_opt, Sol)
end

function solve_qep(qep, x, S, sub_inds; debug=false)
    # TODO why not just preemptively add all subpiece sets to the dictionary,
    # and only modify the qp constraint dependencies? This requires a ton of
    # extra copying of huge matrices for large problems. Too lazy to fix now.
    qep_augmented = deepcopy(qep)
    for qp in values(qep_augmented.qps)
        qp.S[-1] = 1.0 # sid -1 is reserved for this shared set
    end
    qep_augmented.sets[-1] = S

    qp_fair = QP(fair_obj(qep), Dict(id=>sum(qp.S[id] for qp in values(qep_augmented.qps) if id ∈ keys(qp.S); init=0.0) for id in keys(qep_augmented.sets)), sub_inds)
    filter!(p->!iszero(p.second), qp_fair.S)
    qep_augmented.qps[-1] = qp_fair
    solve_qep(qep_augmented, x; debug)
end
