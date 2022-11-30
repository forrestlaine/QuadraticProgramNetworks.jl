@testset "my example 1" begin
    """
    variables := x1 x2 x3 x4
    f1: ½ (x1-x2)²
    f2: ½ || [x1; x2] - [x3; x4] ||²
    """
    Q1 = sparse( [ 1 -1 0 0;
                   -1 1 0 0;
                   0 0 0 0;
                   0 0 0 0.0] )
    Q2 = sparse( [ 1 0 -1 0;
                   0 1 0 -1;
                   -1 0 1 0;
                   0 -1 0 1.0] )
    q1 = q2 = zeros(4)
    f1 = Quadratic(Q1, q1)
    f2 = Quadratic(Q2, q2)

    A = sparse([1.0 0 0 0;])
    l = [0.0,]
    u = [Inf,]
    S = Poly(A, l, u)
    sets = Dict(1=>S)
    qp1 = QP(f1, Dict(1=>1.0), [1,])
    qp2 = QP(f2, Dict{Int, Float64}(), [2,])
    qps = Dict(1=>qp1, 2=>qp2)
    net = [Set([2,]), Set([1,])] 
    qp_net = QPNet(qps, sets, net)
    x = zeros(4)
    x, Sol = solve(qp_net, x; debug=true, gen_Sol=true)
    Sol = collect(Sol)
    @test [0.0; 0; -2; 1] ∈ Sol
    @test [1.0; 1; 1.5; 0.5] ∈ Sol
    @test [0.0; -3; -1; -3] ∈ Sol
    @test [0.0; -3; 4; -3] ∈ Sol
    @test [0.5; 0.5; 4; -3] ∈ Sol
    @test [0.0; 0.0; 0; 0] ∈ Sol
end
