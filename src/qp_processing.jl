function solve_qp(Q, q, A, l, u; tol=1e-8, debug=false, solver=:OSQP)
    if solver == :OSQP
        m = OSQP.Model()
        OSQP.setup!(m; P=sparse(Q), q, A=sparse(A), l, u, verbose=false, polish=true, eps_abs=tol, eps_rel=tol, max_iter=20000)
        ret = OSQP.solve!(m)
        @infiltrate debug
        if ret.info.status_val ∉ (1,2)
            error("Solver failure. Status value is $(ret.info.status_val) ($(ret.info.status)).")
        else
            ret.x
        end
    elseif solver == :PATH
        PATHSolver.c_api_License_SetString("2830898829&Courtesy&&&USR&45321&5_1_2021&1000&PATH&GEN&31_12_2025&0_0_0&6000&0_0")
        n = size(Q,1)
        m = size(A,1)
        M = [Q -A' spzeros(m,m);
             A spzeros(m,m) -sparse(1.0I, m,m)
             spzeros(m,m) sparse(1.0I, m, m) spzeros(m,m)] |> SparseMatrixCSC{Float64, Int32}
        qq = [q; zeros(m); zeros(m)] 
        ll = [fill(-Inf, n+m); l]
        uu = [fill(+Inf, n+m); u]
        (path_status, z, info) =  PATHSolver.solve_mcp(M, qq, ll, uu, zeros(2m+n), 
                                                   silent=true, 
                                                   convergence_tolerance=1e-8, 
                                                   cumulative_iteration_limit=100000,
                                                   restart_limits=5,
                                                   lemke_rank_deficiency_iterations=1000)
        @infiltrate debug
        if path_status != PATHSolver.MCP_Solved
            error("Solver failure. Status value is $(path_status)")
        else
            z[1:n]
        end
    else
        error("Solver not supported")
    end
end

function check_qp_convexity(Q, A, l, u, dec_inds, id; tol=1e-6, debug=false)
    p = Poly(A, l, u)
    (; implicitly_equality, vals) = implicit_bounds(p; tol, debug)
    Ae = collect(A[implicitly_equality,dec_inds])
    F = svd(Ae; full=true)
    V = F.V
    r = rank(Diagonal(F.S))
    Z = V[:,r+1:end]
    QQ = collect(Z'*Q[dec_inds, dec_inds]*Z)
    E = eigen(QQ)
    convex = all(E.values .> -tol) 
    if !convex
        error("QP $id is not convex. Exiting.")
    end
end

function verify_solution(qp, id, constraints, dec_inds, x, check_convexity; tol=1e-4, debug=false)
    Q = qp.f.Q[dec_inds,:]
    q = qp.f.q[dec_inds]
    q̃ = Q*x + q

    A, l, u = isempty(constraints) ? (spzeros(0,length(x)), zeros(0), zeros(0)) :
    map(constraints) do c
        (; A, l, u) = vectorize(c)
        A, l, u
    end |> (x -> vcat.(x...))
    m = size(A,1)

    check_convexity && check_qp_convexity(qp.f.Q, A, l, u, dec_inds, id)
    
    ax = A*x 

    feasible = all(in(x, P; tol=1e-3) for P in constraints)
    #feasible = all(ax .> l.-tol) && all(ax .< u.+tol)
    if !feasible
        return (; solution=false, λ = nothing, e="Current point is infeasible when using tolerance $tol.", debug_data=nothing)
    else
        if m == 0
            if isapprox(q̃, zero(q̃); atol=tol)
                return (; solution=true, λ=zeros(0), debug_data=nothing)
            else
                return (; solution=false, λ=nothing, e="Current point is suboptimal", debug_data=nothing)
            end
        else
            pos_inds = ax .< l .+ 1e-2
            neg_inds = ax .> u .- 1e-2

            both_inds = pos_inds .&& neg_inds
            pos_inds = pos_inds .&& .!both_inds
            neg_inds = neg_inds .&& .!both_inds

            A₊ = A[pos_inds,dec_inds]
            A₋ = A[neg_inds,dec_inds]
            A₀ = A[both_inds,dec_inds]

            n₊ = sum(pos_inds)
            n₋ = sum(neg_inds)

            try

                Ā = [A₊' -A₋' A₀']
                λ = Ā\q̃
                λ₊ = λ[1:n₊]
                λ₋ = λ[n₊+1:n₊+n₋]
                λ₀ = λ[n₊+n₋+1:end]
                if all(λ₊ .> -tol) && all(λ₋ .> -tol) && isapprox(Ā*λ, q̃; atol=tol)
                    λ_out = zeros(m)
                    λ_out[pos_inds] = λ₊
                    λ_out[neg_inds] = -λ₋
                    λ_out[both_inds] = λ₀
                    return (; solution=true, λ = λ_out)
                else
                    error("Simple solve didn't work.")
                end
            catch e
                lb = map(1:m) do i
                    neg_inds[i] || both_inds[i] ? -Inf : 0.0
                end
                ub = map(1:m) do i
                    pos_inds[i] || both_inds[i] ? +Inf : 0.0
                end
                Ad = A[:,dec_inds]
                try
                    λ = solve_qp(Ad*Ad', -Ad*q̃, sparse(I, m, m), lb, ub; solver=:PATH)
                    if isapprox(Ad'*λ, q̃; atol=1e-4)
                        return (; solution=true, λ, debug_data=nothing)
                    else
                        return (; solution=false, λ, e="Current point is suboptimal (via QP).", debug_data=(; Ad, q̃, λ, pos_inds, neg_inds, both_inds, ax, l, u))
                    end
                catch ee
                    return (; solution=false, λ=nothing, e="Solving for duals failed. $ee", debug_data=(; ee))
                end
            end
        end
    end 
end
    
function process_qp(qpn::QPNet, id::Int, x, S; exploration_vertices=0)
    qp = qpn.qps[id]
    base_constraints = [qpn.constraints[c].poly for c in qp.constraint_indices]
    dec_inds = decision_inds(qpn, id)
    check_convexity = qpn.options.check_convexity

    Solgraphs_out = Dict()
    gen_solution_graphs = !(id in qpn.network_depth_map[1]) || qpn.options.gen_solution_map

    child_inds = qpn.network_edges[id] |> collect
    if length(child_inds) > 0
        cardinalities = map(child_inds) do j
            1:length(S[j])
        end
        if any(length.(cardinalities) .< 1) 
            error("Solution graphs were not properly populated.")
        end
        
        all_subpiece_indices = collect(Iterators.product(cardinalities...))
        @debug "Creating subpiece solgraphs in the processing of node $id. 
                Number of combinations of subpieces to investigate: $(length(all_subpiece_indices))."
        subpiece_solgraphs = map(all_subpiece_indices) do child_subpiece_indices
            let child_inds=child_inds, 
                child_subpiece_indices=child_subpiece_indices, 
                S=S, 
                base_constraints=base_constraints, 
                qp=qp, 
                dec_inds=dec_inds, 
                x=x, 
                exploration_vertices=exploration_vertices, 
                check_convexity=check_convexity

                children_solgraph_polys = map(child_inds, child_subpiece_indices) do j,ji
                    P = S[j][ji]
                end
                appended_constraints = [base_constraints; children_solgraph_polys]
                ret = verify_solution(qp, id, appended_constraints, dec_inds, x, check_convexity)
                if !ret.solution
                    subpiece_assignments = Dict(j=>ji for (j,ji) in zip(child_inds, child_subpiece_indices))
                    (; solution=false, e=ret.e, subpiece_assignments, ret.debug_data)
                else
                    if gen_solution_graphs
                        solgraph_generator = process_solution_graph(qp, appended_constraints, dec_inds, x, ret.λ; exploration_vertices)
                        high_dim = length(solgraph_generator.z) + length(solgraph_generator.w)
                        @debug "There are $(length(solgraph_generator.unexplored_Ks)) nodes to expand, excluding additional vertex exploration. 
                                Projecting from $high_dim to $(length(x)) dimensions. Exploration vertices is $exploration_vertices"
                        solgraph = (children_solgraph_polys, remove_subsets(PolyUnion(collect(solgraph_generator))))
                        solgraph
                    else
                        solgraph = nothing
                    end
                    (; solution=true, solgraph)
                end
            end 
        end
        results = fetch.(subpiece_solgraphs)
        for r in results
            if !r.solution
                @debug "When using one of the subpieces, the current point is not a solution for QP $id."
                @debug r.e
                return (; solution=false, r.e, failed=false, r.subpiece_assignments, r.debug_data)
            end
        end
        @debug "The current point is optimal for QP $id using each of the subpieces. 
                Solgraphs have been found for each subpiece. Going to now combine them."
        if gen_solution_graphs
            try
                S_out = combine((r.solgraph for r in results), x; show_progress=true) |> collect |> PolyUnion
            catch e
                return (; solution = false, failed=true, S=nothing)
            end
        else
            S_out = nothing
        end
    else
        ret = verify_solution(qp, id, base_constraints, dec_inds, x, check_convexity)
        if !ret.solution
            return (; solution=false, ret.e, failed=false, subpiece_assignments = Dict())
        else
            if gen_solution_graphs
                S_out = process_solution_graph(qp, base_constraints, dec_inds, x, ret.λ; exploration_vertices) |> collect |> PolyUnion
                if length(S_out) == 0
                    error("This shouldn't happen. Solution graph is empty.")
                end
            else
                S_out = nothing
            end
        end
    end
    return (; solution=true, S=S_out, failed=false)
end 

function combine(solgraphs, x; show_progress=true)
    regions = Poly[]
    solutions = PolyUnion[]
    for (r,s) in solgraphs
        local pr
        try
            pr = poly_intersect(r...)
        catch e
            @infiltrate
        end
        pr = project(pr, 1:embedded_dim(pr))
        push!(regions,pr)
        push!(solutions,s)
    end
    combine(regions, solutions, x; show_progress)
end

"""
Conustructs the solution set
S := ⋃ₚ ⋂ᵢ Zᵢᵖ

Zᵢᵖ ∈ { Rᵢ', Sᵢ }
where Rᵢ' is the set complement of Rᵢ.
"""
function combine(regions, solutions, x; show_progress=true)
    if length(solutions) == 0
        error("No solutions to combine... length solutions: 0, length regions: $(length(regions))")
    elseif length(solutions) == 1
        PolyUnion(collect(first(solutions)))
    else
        complements = map(complement, regions)
        it = 0
        combined = map(zip(solutions, complements)) do (s, rc)
            it += 1
            PolyUnion([collect(s); rc.polys]) #|> remove_subsets
        end
        widths = [length(c) for c in combined]

        if length(widths) > 3 && sum(widths) > 20
            @error("Too many solutions to combine. $widths")
            error("Too many solutions to combine.")
        end
        @debug "Widths: $widths"
        #combined = [[collect(s); rc] for (s, rc) in zip(solutions, complements)]
        root = IntersectionRoot(combined, length.(complements), x; show_progress)
        root
    end
end
