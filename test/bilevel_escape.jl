@testset "simple bilevel" begin
    qpn = setup(:bilevel_escape; make_requests=true)
    ret = solve(qpn, zeros(4))
    @test ret.x_opt ≈ [2.0, 0, 1, 0]

    ret = solve(qpn, [-1.0, -1, 0, 0])
    @test ret.x_opt ≈ [2.0, 0, 1, 0]
end
