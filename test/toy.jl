@testset "my example 1" begin
    """
    variables := x1 x2 x3 x4
    f1: ½ (x1-x2)²
    f2: ½ || [x1; x2] - [x3; x4] ||²
    """
    high_dim = true
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
    if high_dim
        xinit = [0.0, 0, -2, 1]
        x, Sol = solve(qp_net, xinit; debug=false, gen_Sol=true, high_dim)
        @test x ≈ [0.0, 0, -2, 1]

        xinit = [0, 0, 1.5, 0.5]
        x, Sol = solve(qp_net, xinit; debug=false, gen_Sol=true, high_dim)
        @test x ≈ [1, 1, 1.5, 0.5]

        xinit = [0.0, 0, -1, -3]
        x, Sol = solve(qp_net, xinit; debug=false, gen_Sol=true, high_dim)
        @test x ≈ [0.0, -3, -1, -3]
        
        xinit = [0.0, 0, 4, -3]
        x, Sol = solve(qp_net, xinit; debug=false, gen_Sol=true, high_dim)
        @test (x ≈ [0.5, 0.5, 4, -3]) || (x ≈ [0.0, -3, 4, -3])
        
        xinit = [0.0, 0, 0, 0]
        x, Sol = solve(qp_net, xinit; debug=false, gen_Sol=true, high_dim)
        @test x ≈ [0.0, 0,0,0]
    else
        xinit = [0.0, 0, 0, 0]
        x, Sol = solve(qp_net, xinit; debug=false, gen_Sol=true, high_dim)
        Sol = collect(Sol)
        @test [0.0; 0; -2; 1] ∈ Sol
        @test [1.0; 1; 1.5; 0.5] ∈ Sol
        @test [0.0; -3; -1; -3] ∈ Sol
        @test [0.0; -3; 4; -3] ∈ Sol
        @test [0.5; 0.5; 4; -3] ∈ Sol
        @test [0.0; 0.0; 0; 0] ∈ Sol
        @test [1.0; 0.5; 4; -3] ∉ Sol
    end
end
