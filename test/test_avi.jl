@testset "LCP examples" begin
    M1 = [0 1 -1;
          -1 0 1;
          1 -1 0]
    M2 = -M1
    e = ones(3)
    Z = zeros(3,3)
    M = [Z M1 -e zeros(3);
         M2' Z zeros(3) -e;
         e' zeros(1,5)
         zeros(1,3) e' zeros(1,2)]
    o = [zeros(6); -1; -1]
    N = zeros(8,0)
    w = zeros(0)
    l = [zeros(6); fill(-Inf, 2)]
    u = [fill(Inf, 6); fill(Inf,2)]
    avi = AVI(sparse(M), sparse(N), o, l, u)
    z0 = rand(8)
    ret = solve_avi(avi, z0, zeros(0))
    @test ret.status == QPN.SUCCESS
end
