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
                                                   silent=false, 
                                                   convergence_tolerance=1e-8, 
                                                   cumulative_iteration_limit=100000,
                                                   restart_limits=5,
                                                   lemke_rank_deficiency_iterations=1000)
    (; sol_bad, degree, r) = check_avi_solution(avi, z, w)
    if sol_bad
        @infiltrate
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
        @infiltrate
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

function relax_gavi(gavi::GAVI, relaxable_inds)
    param_inds = setdiff(1:size(gavi.N,2), relaxable_inds)
    d1 = length(gavi.l1)
    d2 = length(gavi.l2)
    dr = length(relaxable_inds)
    M = [spzeros(dr, d1+d2+dr);
         gavi.N[:,relaxable_inds] gavi.M]
    N = [spzeros(dr, length(param_inds)); gavi.N[:, param_inds]]
    o = [zeros(dr); gavi.o]
    l1 = [fill(-Inf, dr); gavi.l1]
    u1 = [fill(Inf, dr); gavi.u1]

    A = [gavi.B[:,relaxable_inds] gavi.A]
    B = gavi.B[:,param_inds]

    GAVI(M,N,o,l1,u1,A,B,gavi.l2,gavi.u2)
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

function create_labeled_gavi_from_qp(qp_net, id, solution_graphs)
    dvars = decision_inds(qp_net, id) 
    n = length(dvars)
    qp = qp_net.qps[id]
    n_total = size(qp.f.Q,2)
    
    # Z = [x; ξᵢ; λᵢ; ψᵢ]
    # cons = [M1*Z + q1 ⟂ l1 ≤ ξᵢ   ≤ u1;
    #        [[λᵢ; ψᵢ]  ⟂ l2 ≤ M2*Z ≤ u2]
    
    labels = Dict{String, Int}()
    for i in 1:n_total
        labels["x_$i"] = i
    end
    for (e,i) in enumerate(dvars)
        labels["ξ_$(id)_$(i)"] = n_total+e
    end
    total = n_total+n


    #A, l, u = mapreduce(vcat, qp.constraint_indices) do ci
    #    (; A, l, u) = vectorize(qp_net.constraints[ci].poly)
    #    for i in 1:size(A,1)
    #        labels["λ_$(id)_$(ci)_$i"] = total+i
    #    end
    #    total += size(A,1)
    #    A, l, u
    #end |> (x->vcat.(x...))
    
    A,l,u = map(qp.constraint_indices) do ci
        (; A, l, u) = vectorize(qp_net.constraints[ci].poly)
        for i in 1:size(A,1)
            labels["λ_$(id)_$(ci)_$i"] = total+i
        end
        total += size(A,1)
        A, l, u
    end |> (x->vcat.(x...))

    A_Si, l_Si, u_Si = isempty(qp_net.network_edges[id]) ? (spzeros(0, n_total), zeros(0), zeros(0)) :
    map(collect(qp_net.network_edges[id])) do j
        (; A, l, u) = vectorize(solution_graphs[j])        
        for i in 1:size(A,1)
            labels["ψ_$(id)_$(j)_$i"] = total+i
        end
        total += size(A,1)
        A, l, u
    end |> (x->vcat.(x...))
    
    M1 = [qp.f.Q[dvars, :] -sparse(I, n, n) -A[:,dvars]' -A_Si[:,dvars]']
    q1 = qp.f.q[dvars]
    M2 = [A; A_Si]
    l2 = [l; l_Si]
    u2 = [u; u_Si]

    (; dvars, labels, M1, q1, M2, l2, u2)
end


"""
Generalized AVI setup for Nash Equilibrium Problem.

Z := [dvars; [ξᵢ for i in pool];  [[λᵢ; ψᵢ] for i in pool]]

"""
function combine_gavis(n, dec_inds, param_inds, labeled_gavis)

    nd = length(dec_inds)
    total_dual_dim = 0
    total_ξ_dim = 0
    for (id, lgavi) in labeled_gavis
        total_dual_dim += (size(lgavi.M1, 2) - n)
        total_ξ_dim += (size(lgavi.M1,1))
    end
    ξ_offset_ranges = Dict{Int, AbstractVector{Int}}()
    λψ_offset_ranges = Dict{Int, AbstractVector{Int}}()
    offset1 = 0
    offset2 = total_ξ_dim

    player_pool = collect(keys(labeled_gavis))

    M, N, q = map(player_pool) do id

        lgavi = labeled_gavis[id]

        M = lgavi.M1
        q = lgavi.q1

        dual_dim = size(M,2) - n
        ξ_dim = size(M,1)
        λψ_dim = dual_dim - ξ_dim

        ξ_offset_ranges[id] = offset1+1:offset1+ξ_dim
        λψ_offset_ranges[id] = offset2+1:offset2+λψ_dim   

        Mi = spzeros(size(M,1), nd+total_dual_dim)
        Mi[:,1:nd] = M[:, dec_inds]

        Mi[:,nd.+ξ_offset_ranges[id]] = M[:,n+1:n+ξ_dim]
        Mi[:,nd.+λψ_offset_ranges[id]] = M[:,n+ξ_dim+1:end]
        Ni = M[:,param_inds]

        offset1 += ξ_dim
        offset2 += λψ_dim

        Mi, Ni, q  
    end |> (x->vcat.(x...))
    
    A, B, l2, u2 = map(player_pool) do id
        lgavi = labeled_gavis[id]
        A = lgavi.M2
        l = lgavi.l2
        u = lgavi.u2
        A[:,dec_inds], A[:,param_inds], l, u
    end |> (x->vcat.(x...))

    top_M = spzeros(nd, size(M,2))
    top_N = spzeros(nd, size(N,2))
    top_q = zeros(nd)

    for (id, lgavi) in labeled_gavis
        ξ_offset_range = ξ_offset_ranges[id]
        for (di, d) in enumerate(lgavi.dvars)
            top_M[di, nd+ξ_offset_range[lgavi.labels["ξ_$(id)_$d"]-n]] = 1.0
        end
    end

    M = [top_M; M]
    N = [top_N; N]
    o = [top_q; q]
    l1 = fill(-Inf, length(o))
    u1 = fill(+Inf, length(o))
    A = [A spzeros(size(A,1), total_dual_dim)]

    gavi = GAVI(M, N, o, l1, u1, A, B, l2, u2)
end

"""
Solve the Quadratic Equilibrium Problem.
"""
function solve_qep(qp_net, level, x, request, relaxable_inds, S=Dict{Int, Poly}(), shared_decision_inds=Vector{Int}();
                   var_indices=nothing,
                   subpiece_index=0,
                   debug=false,
                   request_comes_from_parent=false,
                   high_dimension=false,
                   make_requests=false,
                   gen_sol=true,
                   shared_variable_mode=SHARED_DUAL,
                   rng=MersenneTwister(1))

    x_dim = length(x)
    player_pool = qp_net.network_depth_map[level]

    dec_inds = union((decision_inds(qp_net, id) for id in player_pool)...) |> sort
    param_inds = setdiff(1:x_dim, dec_inds)

    labeled_gavis = Dict(id=>create_labeled_gavi_from_qp(qp_net, id, S) for id in player_pool)
    gavi = combine_gavis(x_dim, dec_inds, param_inds, labeled_gavis)

    w = x[param_inds]

    #       decision_vars       aux_vars                ψ vars                      dual_vars
    z0 = [x[dec_inds]; zeros(size(gavi.M,2)-length(dec_inds))]

    # TODO : if I want to support successive minimization of ||ψ|| over
    # iterations, need to properly warmstart with previous solution? Might
    # require reducing dimension AFTER psi minimization

    (; z, status) = solve_gavi(gavi, z0, w)

    @infiltrate # Code runs up until this new point. Need to change below

    if status != SUCCESS
        relaxable_parent_inds = setdiff(relaxable_inds, decision_inds)
        relaxable_parent_inds = [findfirst(param_inds .== i) for i in relaxable_parent_inds]
        if !isempty(relaxable_parent_inds) && request_comes_from_parent
            debug && @info "AVI solve error, but some parameter variables can be relaxed. Constructing relaxed GAVI."
            r_gavi = relax_gavi(gavi, relaxable_parent_inds)
            r_z0 = [w[relaxable_parent_inds]; z0]
            r_w = w[setdiff(1:length(w), relaxable_parent_inds)]
            ret = solve_gavi(r_gavi, r_z0, r_w)
            status = ret.status
            if status != SUCCESS
                error("AVI solve error, even after relaxing indices.")
            end
            l_r = length(relaxable_parent_inds)
            w[relaxable_parent_inds] = ret.z[1:l_r]
            z = ret.z[l_r+1:end]
        else
            @infiltrate
            error("AVI solve error!")
        end
    end
    

    ψ_inds = collect(N_private_vars+N_shared_vars+1:N_private_vars+N_shared_vars*(N_players+1))

    if high_dimension
        throw(error("High dimension mode not supported at the moment"))
        extra_rounds = level==1 ? 0 : 5
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
            x_opt[param_inds] = w

            # TODO figure out request structure with vertex expansion (is
            # v-enum even required?)
            Sol = gen_sol ? LocalGAVISolutions(gavi, z, w, level, subpiece_index, decision_inds, param_inds, request; max_vertices = 1000) : nothing
            @debug "Solution map generated."

            # TODO : should probably propagate any parent level requests if
            # they appear in S 

            if isnothing(S) || !make_requests
                identified_request = Set{Linear}()
            else
                S_duals = z[S_dual_inds]
                identified_request = identify_request(S, S_duals, request; propagate=request_comes_from_parent)
            end
            (; x_opt, Sol, identified_request)
        else
            @error "Invalid shared variable mode: $shared_variable_mode."
        end
    end
end

function identify_request(S, λs, parent_request; propagate=false)
    identified_request = Set{Linear}()
    (; A, l, u) = vectorize(S)
    (m,d) = size(A)

    if propagate
        for req in parent_request
            if iszero(req.a[d+1:end])
                for i = 1:m
                    if req.a[1:d] ≈ A[i,:] 
                        union!(identified_request, propagate_request(A[i,:], get_parent(S, i)))
                    elseif req.a[1:d] ≈ -A[i,:]
                        union!(identified_request, propagate_request(-A[i,:], get_parent(S, i)))
                    end
                end
            end
        end
    else
        for (i, λ) in enumerate(λs)
            if λ ≥ 1e-4 && has_parent(S, i)
                union!(identified_request, propagate_request(A[i,:], get_parent(S, i)))
            elseif λ ≤ -1e-4 && has_parent(S, i)
                union!(identified_request, propagate_request(-A[i,:], get_parent(S, i)))
            end
        end
    end
    identified_request
end

function propagate_request(request, poly)
    m = OSQP.Model()
    d = embedded_dim(poly)
    n = length(request)
    q = zeros(d)
    q[1:n] = request
    (; A, l, u) = vectorize(poly)
    OSQP.setup!(m; q, A, l, u,
                verbose=false,
                polish=true,
                eps_abs=1e-8,
                eps_rel=1e-8)
    ret = OSQP.solve!(m)
    prop_requests = Set{Linear}()
    if ret.info.status_val == 1
        duals = -ret.y
        for (i, λ) in enumerate(duals)
            if λ ≥ 1e-4
                push!(prop_requests, Linear(A[i,:]))
            elseif λ ≤ -1e-4
                push!(prop_requests, Linear(-A[i,:]))
            end
        end
    else
        # This shouldn't happen (would mean halfspace in projected poly isn't
        # implied by halfspace in parent poly)
        throw(error("Unable to propagate request to parent poly for some reason."))
    end
    prop_requests
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
