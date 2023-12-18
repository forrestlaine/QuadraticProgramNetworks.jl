function solve_qp(Q, q, A, l, u; tol=1e-8)
    m = OSQP.Model()
    OSQP.setup!(m; P=sparse(Q), q, A=sparse(A), l, u, verbose=false, polish=true, eps_abs=tol, eps_rel=tol)
    ret = OSQP.solve!(m)
    if ret.info.status_val != 1
        error("Solver failure. Status value is $(ret.info.status_val).")
    else
        ret.x
    end
end

function verify_solution(qp, constraints, dec_inds, x; tol=1e-6)
    Q = qp.f.Q[dec_inds,:]
    q = qp.f.q[dec_inds]
    q̃ = Q*x + q

    A, l, u = isempty(constraints) ? (spzeros(0,length(x)), zeros(0), zeros(0)) :
    map(constraints) do c
        (; A, l, u) = vectorize(c)
        A, l, u
    end |> (x -> vcat.(x...))
    m = size(A,1)
    
    ax = A*x 

    feasible = all(ax .> l.-tol) && all(ax .< u.+tol)
    if !feasible
        return (; solution=false, λ = nothing, e="Current point is infeasible when using tolerance $tol.")
    else
        if m == 0
            if isapprox(q̃, zero(q̃); atol=tol)
                return (; solution=true, λ=zeros(0))
            else
                return (; solution=false, λ=nothing, e="Current point is suboptimal")
            end
        else
            pos_inds = ax .< l .+ tol
            neg_inds = ax .> u .- tol
            both_inds = pos_inds .&& neg_inds
            pos_inds = pos_inds .&& .!both_inds
            neg_inds = neg_inds .&& .!both_inds

            A₊ = A[pos_inds,dec_inds]
            A₋ = A[neg_inds,dec_inds]
            A₀ = A[both_inds,dec_inds]

            n₊ = sum(pos_inds)
            n₋ = sum(neg_inds)

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
                try 
                    lb = map(1:m) do i
                        neg_inds[i] || both_inds[i] ? -Inf : 0.0
                    end
                    ub = map(1:m) do i
                        pos_inds[i] || both_inds[i] ? +Inf : 0.0
                    end
                    Ad = A[:,dec_inds]
                    λ = solve_qp(Ad*Ad', -Ad*q̃, sparse(I, m, m), lb, ub)
                    if isapprox(A'*λ, q̃; atol=tol)
                        return (; solution=true, λ)
                    else
                        return (; solution=false, λ=nothing, e="Current point is suboptimal")
                    end
                catch e
                    @infiltrate
                    return (; solution=false, λ=nothing, e)
                end
            end
        end
    end 
end
    

function process_qp(qpn::QPNet, id::Int, x)
    qp = qpn.qps[id]
    constraints = [qpn.constraints[c].poly for c in qp.constraint_indices]
    dec_inds = decision_inds(qpn, id)
    
    ret = verify_solution(qp, constraints, dec_inds, x)
    if !ret.solution
        return (; solution=false, ret.e)
    else
        S = process_solution_graph(qp, constraints, dec_inds, x, ret.λ) |> collect |> PolyUnion
        return (; solution=true, S)
    end
end

function process_qp(qpn::QPNet, id::Int, x, S)
    qp = qpn.qps[id]
    base_constraints = [qpn.constraints[c].poly for c in qp.constraint_indices]
    dec_inds = decision_inds(qpn, id)

    Solgraphs_out = Dict()

    child_inds = qpn.network_edges[id] |> collect
    if length(child_inds) > 0
        cardinalities = map(child_inds) do j
            1:length(S[j])
        end
        for child_subpiece_indices in Iterators.product(cardinalities...)
            children_solgraph_polys = map(child_inds, child_subpiece_indices) do j,ji
                P = S[j][ji]
            end
            appended_constraints = [base_constraints; children_solgraph_polys]
            ret = verify_solution(qp, appended_constraints, dec_inds, x)
            if !ret.solution
                subpiece_assignments = Dict(j=>ji for (j,ji) in zip(child_inds, child_subpiece_indices))
                return (; solution=false, subpiece_assignments)
            else
                Solgraphs_out[child_subpiece_indices] = (children_solgraph_polys, PolyUnion(collect(process_solution_graph(qp, appended_constraints, dec_inds, x, ret.λ))))
            end
        end
        S_out = combine(values(Solgraphs_out), x) |> collect |> PolyUnion # Maybe this can stay an iterator, but for now I am going to collect

    else
        ret = verify_solution(qp, base_constraints, dec_inds, x)
        if !ret.solution
            return (; solution=false, ret.e)
        else
            S_out = process_solution_graph(qp, base_constraints, dec_inds, x, ret.λ) |> collect |> PolyUnion
        end
    end

    return (; solution=true, S=S_out)
end 

function combine(solgraphs, x; show_progress=true)
    regions = Poly[]
    solutions = PolyUnion[]
    for (r,s) in solgraphs
        pr = poly_intersect(r...)
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
            PolyUnion([collect(s); rc.polys])
        end
        #combined = [[collect(s); rc] for (s, rc) in zip(solutions, complements)]
        IntersectionRoot(combined, length.(complements), x; show_progress)
    end
end
