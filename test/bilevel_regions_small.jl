@testset "bilevel regions toy" begin
    """ 
    QP1: min 0.5(x-w)^2
          x
    QP2: max x*y
          0 ≤ y ≤ 1

    variables w x y
    """
    Q1 = sparse( [ 1 -1 0;
                   -1 1 0;
                   0 0 0.0] )
    Q2 = sparse( [ 0 0 0;
                   0 0 -1;
                   0 -1 0.0] )
    q1 = q2 = zeros(3)
    f1 = Quadratic(Q1, q1)
    f2 = Quadratic(Q2, q2)

    A = sparse([0 0 1.0])
    l = [0.0,]
    u = [1.0,]
    S = Poly(A, l, u)

    sets = Dict(1=>S)

    qp1 = QP(f1, Dict{Int, Float64}(), [2,])
    qp2 = QP(f2, Dict(1=>1.0), [3,])
    qps = Dict(1=>qp1, 2=>qp2)
    net = [Set([1,]), Set([2,])] 
    qp_net = QPNet(qps, sets, net)

    x = [1.0, -1.0, 1.0] 
    x, Sol = solve(qp_net, x; debug=true)
    @test [1.0, 1.0, 1.0] ≈ x 
end
