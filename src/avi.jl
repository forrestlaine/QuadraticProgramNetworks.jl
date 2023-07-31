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
    (path_status, z, info) =  PATHSolver.solve_mcp(avi.M, avi.N*w+avi.o,avi.l, avi.u, z0, 
                                                   silent=true, 
                                                   convergence_tolerance=1e-8, 
                                                   cumulative_iteration_limit=100000,
                                                   restart_limits=5,
                                                   lemke_rank_deficiency_iterations=1000)
    (; sol_bad, degree, r) = check_avi_solution(avi, z, w)
    if sol_bad
        return (; z, status=FAILURE)
    end
    status = (path_status == PATHSolver.MCP_Solved || path_status == PATHSolver.MCP_Solved) ? SUCCESS : FAILURE
    return (; z, status)
end

function find_closest_feasible!(gavi, z0, w)
    model = OSQP.Model()
    n = length(z0)
    c = gavi.B*w
    OSQP.setup!(model;
                P = sparse(I, n, n),
                q = -z0,
                A = gavi.A,
                l = gavi.l2-c,
                u = gavi.u2-c,
                verbose=false,
                polish=true,
                eps_abs=1e-8,
                eps_rel=1e-8)
    ret = OSQP.solve!(model)
    if ret.info.status_val == 1
        z0 .= ret.x
    else
        @warn "Feasible initialization not cleanly solved. Solve status: $(ret.info.status)"
    end
end

function solve_gavi(gavi::GAVI, z0, w; presolve=true)
    presolve && find_closest_feasible!(gavi, z0, w)
    avi = convert(gavi)
    d1 = length(gavi.l1)
    d2 = length(gavi.l2)
    s = gavi.A*z0+gavi.B*w
    z0s = copy([z0; s])
    (; z, status) = solve_avi(avi, z0s, w)
    zg = z[1:d1+d2]
    (; z=zg, status)
end

function convert(gavi::GAVI)
    d1 = length(gavi.l1)
    d2 = length(gavi.l2)
    M = [gavi.M spzeros(d1,d2);
         gavi.A -sparse(I,d2,d2);
         spzeros(d2,d1) sparse(I, d2,d2) spzeros(d2,d2)]
    # [M   0]          z1
    # [A  -I] z  ⟂     z2
    # [0I  0]         l≤s≤u 

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
    bad_count = sum(abs.(z[r_pos]-avi.l[r_pos]) .> tol) + 
                sum(abs.(z[r_neg]-avi.u[r_neg]) .> tol) +
                sum(z-avi.l .< -tol) + sum(z-avi.u .> tol)
    return (; sol_bad = bad_count > 0, degree = bad_count, r)
end

"""
Solve the Quadratic Equilibrium Problem.
"""
function solve_qep(qep_base, x, request, S=nothing, shared_decision_inds=Vector{Int}();
                   var_indices=nothing,
                   level=0,
                   subpiece_index=0,
                   debug=false,
                   high_dimension=false,
                   gen_sol=true,
                   shared_variable_mode=SHARED_DUAL,
                   rng=MersenneTwister(1))

    x_dim = length(x)
    N_players = length(qep_base.qps)

    if isnothing(S)
        qep = qep_base
        aux_dim = 0
        underconstrained=false
        @assert length(shared_decision_inds) == 0
    else
        qep = deepcopy(qep_base)
        foreach(qp->push!(qp.constraint_indices, -1), values(qep.qps))
        if shared_variable_mode==SHARED_DUAL
            qep.constraints[-1] = Constraint(S, Dict(i=>1 for i in keys(qep.qps))) 
        elseif shared_variable_mode==MIN_NORM
            qep.constraints[-1] = Constraint(S, Dict(i=>i for i in keys(qep.qps))) 
        end
        aux_dim = embedded_dim(S) - x_dim
        constrained_lower = embedded_dim(S) - intrinsic_dim(S)
        underconstrained = constrained_lower < (aux_dim + length(shared_decision_inds))
        @assert length(shared_decision_inds) > 0
    end
    player_order = collect(keys(qep.qps)) |> sort
    constraint_order = collect(keys(qep.constraints)) |> sort

    private_decision_inds = reduce(vcat, (qep.qps[k].var_indices for k in player_order)) #|> sort
    decision_inds = [private_decision_inds; shared_decision_inds]
    
    N_shared_vars = length(shared_decision_inds) + aux_dim
    N_private_vars = sum(length.(private_decision_inds))
   
    standard_dual_dim = sum(length(Set(values(C.group_mapping)))*length(C.poly) for C in values(qep.constraints))

    param_inds = setdiff(Set(1:x_dim), Set(decision_inds)) |> collect #|> sort
    
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
    

    offset = N_private_vars + (N_players+1)*N_shared_vars
    S_dual_inds = nothing
    As, A2s, Bs, ls, us = map(constraint_order) do id
        # TODO make sure that As / A2s accounts for auxiliary variables
        # properly (should be part of shared vars I think)
        (; A,l,u) = vectorize(qep.constraints[id].poly)
        con_dim = length(l)
        if id == -1
            S_dual_inds = collect(offset+1:offset+con_dim)
        else
            offset += con_dim
        end
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
        group_labels = keys(group_to_player_map) |> collect |> sort
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
            end
        end
        [A1 Ax], A2, B1, l1, u1
    end |> (x->vcat.(x...))
    
    M11 = Qs
    M12 = underconstrained ? Mψ : spzeros(size(Mψ))
    M13 = -A2s'
    M21 = spzeros(N_shared_vars, N_private_vars+N_shared_vars)
    M22 = repeat(sparse(I, N_shared_vars, N_shared_vars), 1, N_players)
    M23 = spzeros(N_shared_vars, standard_dual_dim)
    M31 = As
    M32 = spzeros(standard_dual_dim, N_shared_vars*N_players)
    M33 = spzeros(standard_dual_dim, standard_dual_dim)

    M = [M11 M12 M13;
         M21 M22 M23]

    N1 = Rs
    N2 = spzeros(N_shared_vars, length(param_inds))
    N3 = Bs

    N = [N1; N2]

    o = [qs; zeros(N_shared_vars)]

    l1 = fill(-Inf, length(qs) + N_shared_vars)
    u1 = fill(Inf, length(qs) + N_shared_vars)

    A = [M31 M32 M33]
    B = N3

    l2 = ls
    u2 = us

    w = x[param_inds]
   
    #       decision_vars       aux_vars                ψ vars                      dual_vars
    z0 = [x[decision_inds]; zeros(aux_dim); zeros(N_players*N_shared_vars); zeros(standard_dual_dim)]
    # TODO : if I want to support successive minimization of ||ψ|| over
    # iterations, need to properly warmstart with previous solution? Might
    # require reducing dimension AFTER psi minimization

 
    gavi = GAVI(M,N,o,l1,u1,A,B,l2,u2)
    (; z, status) = solve_gavi(gavi, z0, w)


    #status != SUCCESS && @infiltrate
    status != SUCCESS && @warn "AVI solve error!"
    status != SUCCESS && error("AVI solve error!")
    ψ_inds = collect(N_private_vars+N_shared_vars+1:N_private_vars+N_shared_vars*(N_players+1))

    if high_dimension
        extra_rounds = level==1 ? 0 : 5
        #level == 3 && @infiltrate
        z_orig = z
        (; piece, x_opt, reduced_inds, z) = get_single_solution(gavi,z,w,level,subpiece_index,decision_inds,param_inds,rng; debug=false, permute=false, extra_rounds, level)
        z_inds_remaining = setdiff(1:length(z), reduced_inds)
        z = z[z_inds_remaining] 
        if length(ψ_inds) > 0 && underconstrained && shared_variable_mode == MIN_NORM
            ψ_inds_remaining = setdiff(ψ_inds, reduced_inds)
            f_min_norm = min_norm_objective(length(z), ψ_inds_remaining)
            (; piece, x_opt, z_revised) = revise_avi_solution(f_min_norm, piece, z, w, decision_inds, param_inds, rng)
        end
	    @infiltrate !([z;w] in piece)
	    xz_permuted = zeros(length(z)+length(w))
	    xz_permuted[decision_inds] = z[1:length(decision_inds)]
	    xz_permuted[param_inds] = w
	    xz_permuted[length(decision_inds)+length(param_inds)+1:end] = z[length(decision_inds)+1:end]
        permute!(piece, decision_inds, param_inds)
	    @infiltrate !(xz_permuted in piece)
        reduced_piece = eliminate_variables(piece, x_dim+1:embedded_dim(piece), xz_permuted)
	    @infiltrate !(x_opt in reduced_piece)
        @infiltrate embedded_dim(reduced_piece) > length(x_opt)
        (; x_opt, Sol=[reduced_piece,], identified_request=nothing)
    else
        if shared_variable_mode == MIN_NORM
            @error "not implemented yet" 
        elseif shared_variable_mode == SHARED_DUAL
            @debug "Found solution, now generating solution map (level $(level))"
            x_opt = copy(x)
            x_opt[decision_inds] = z[1:length(decision_inds)]

            # TODO figure out request structure with vertex expansion (is
            # v-enum even required?)
            
            Sol = gen_sol ? LocalGAVISolutions(gavi, z, w, level, subpiece_index, decision_inds, param_inds, request; max_vertices = 0) : nothing
            @debug "Solution map generated."
            if isnothing(S)
                identified_request = Set{Linear}()
            else
                S_duals = z[S_dual_inds]
                @assert length(S_duals) == length(S)
                (; A, l, u) = vectorize(S)
                identified_request = Linear[]
                for (i,λ) in enumerate(S_duals)
                    if λ ≥ 1e-4
                        push!(identified_request, Linear(A[i,:]))
                    elseif λ ≤ -1e-4
                        push!(identified_request, Linear(-A[i,:]))
                    end
                end
                identified_request = Set(identified_request)
            end
            (; x_opt, Sol, identified_request)
        else
            @error "Invalid shared variable mode: $shared_variable_mode."
        end
    end
end

"""
f(z) = 0.5*∑(zᵢ²; i∈inds)
"""
function min_norm_objective(n, inds)
    Q = spzeros(n,n)
    foreach(i->Q[i,i]=1.0, inds) 
    Quadratic(Q, zeros(n))
end

function revise_avi_solution(f, piece, zr, w, decision_inds, param_inds, rng)
    # TODO refactor this to use solve_qep (need to call this function from
    # algorithm.jl)

    vec = vectorize(piece)
    A = vec.A
    ll = vec.l
    uu = vec.u

    (m,n) = size(A)

    nz = length(zr)
    nw = length(w)

    B = A[:,nz+1:nz+nw]
    A = A[:,1:nz]

    M = [f.Q -A']
    N = spzeros(nz,nw)
    o = f.q
    l1 = fill(-Inf, nz)
    u1 = fill(Inf, nz)
    A2 = [A spzeros(m,m)]
    l2 = ll
    u2 = uu
    
    gavi = GAVI(M, N, o, l1, u1, A2, B, l2, u2)
    z0 = [zr; zeros(m)]
    local z, status
    try 
        (; z, status) = solve_gavi(gavi, z0, w)
    catch e
        @infiltrate
    end
    status != SUCCESS && @infiltrate
    status != SUCCESS && error("AVI solve error!")
    (; piece, x_opt, reduced_inds) = get_single_solution(gavi, z, w, level, subpiece_index, decision_inds, param_inds, rng; permute=false)
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
