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

    A, l, u = map(constraints) do c
        (; A, l, u) = vectorize(c)
        A, l, u
    end |> (x -> vcat.(x...))
    m = size(A,1)
    
    ax = A*x 

    feasible = all(ax .> l.-tol) && all(ax .< u.+tol)
    if !feasible
        return (; solution=false, λ = nothing, e=error("Current point is infeasible when using tolerance $tol."))
    else
        pos_inds = ax .< l .+ tol
        neg_inds = ax .> u .- tol
        A₊ = A[pos_inds,:]
        A₋ = A[neg_inds,:]

        n₊ = length(pos_inds)
        n₋ = length(neg_inds)

        Ā = [A₊ -A₋]
        λ = Ā\q̃
        if all(λ .> -tol) && isapprox(Ā*λ, q̃; atol=tol)
            λ_out = zeros(m)
            λ_out[pos_inds] = +λ[1:n₊]
            λ_out[neg_inds] = -λ[n₊+1:end]
            return (; solution=true, λ = λ_out)
        else
            try 
                lb = map(1:m) do i
                    i in neg_inds ? -Inf : 0.0
                end
                ub = map(1:m) do i
                    i in pos_inds ? +Inf : 0.0
                end
                λ = solve_qp(A*A', -A*q̃, sparse(I, m, m), lb, ub)
                return (; solution=true, λ)
            catch e
                return (; solution=false, λ=nothing, e)
            end
        end
    end 
end
    

function process_qp(qpn::QPNet, id::Int, x)
    qp = qpn.qps[id]
    constraints = [qpn.constriants[c].poly for c in qp.constraint_indices]
    dec_inds = decision_inds(qpn, id)
    
    ret = verify_solution(qp, constraints, dec_inds, x)
    if !ret.solution
        return (; solution=false)
    else
        S = process_solution_graph(qp, constraints, dec_inds, x, ret.λ)
        return (; solution=true, S)
    end
end

function process_qp(qpn::QPNet, id::Int, x, S)
    qp = qpn.qps[id]
    base_constraints = [qpn.constriants[c].poly for c in qp.constraint_indices]
    dec_inds = decision_inds(qpn, id)

    Solgraphs_out = Dict()

    child_inds = qpn.network_edges[id] |> collect
    cardinalities = map(child_inds) do j
        1:length(S[j])
    end
    for child_subpiece_indices in Iterators.product(cardinalities...)
        children_solgraph_polys = map(child_inds, child_subpiece_indices) do j,ji
            P = S[j][ji]
        end
        appended_constraints = [constraints; children_solgraph_polys]
        ret = verify_solution(qp, appended_constraints, dec_inds, x)
        if !ret.solution
            subpiece_assignments = Dict(j=>ji for (j,ji) in zip(child_inds, child_subpiece_indices))
            return (; solution=false, subpiece_assignments)
        else
            Solgraphs_out[child_subpiece_indices] = (children_solgraph_polys, process_solution_graph(qp, appended_constraints, dec_inds, x, ret.λ))
        end
    end

    S_out = combine(values(Solgraphs_out), x) |> collect # Maybe this can stay an iterator, but for now I am going to collect
    return (; solution=true, S_out)
end 

function combine(solgraphs, x; show_progress=true)
    regions = []
    solutions = []
    foreach((r,s)->push!(regions, r); push!(solutions,s), solgraphs)
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


     
    


    
end
