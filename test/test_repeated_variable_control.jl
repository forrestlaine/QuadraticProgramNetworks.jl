@testset "repeated_variable_control" begin
    include("../examples/repeated_variable_control.jl")
    (; qp_net, Q, q, A, l, u) = setup()
    ret = solve(qp_net, zeros(length(qp_net.variables)))
    
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
    @test isapprox(ret.x_opt[4:5], zeros(2); atol=1e-4)
end

