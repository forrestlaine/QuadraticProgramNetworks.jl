using OSQP
using Random
using SparseArrays

@testset "repeated_variable_control" begin
    (; qpn, Q, q, A, l, u) = setup(:repeated_variable_control; make_requests=true)
    ret = solve(qpn, zeros(length(qpn.variables)))
    
    mod = OSQP.Model()
    OSQP.setup!(mod;
                q,
                P = sparse(Q),
                A = sparse(A),
                l,
                u,
                polish=true,
                verbose=false)
    qp_ret = OSQP.solve!(mod)

    @test ret.x_opt[1:3] ≈ qp_ret.x
    @test isapprox(ret.x_opt[4], 0.0; atol=1e-4)
end

