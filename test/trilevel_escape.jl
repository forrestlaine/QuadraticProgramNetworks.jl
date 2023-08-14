@testset "trilevel escape" begin
    qpn = setup(Val(:trilevel_escape))
    initializations = [ [-2.0, -0.5, -1, -.5], fill(-0.5, 4), fill(0.0, 4), fill(0.5, 4), [1.0, 1.0, 1.0, 0.5], [2.0, 1.0, 1.0, 0.5]]
    for init in initializations
        ret = solve(qpn, init)
        @test ret.x_opt ≈ [2.0, 1.0, 1.0, 0.5]
    end
end
