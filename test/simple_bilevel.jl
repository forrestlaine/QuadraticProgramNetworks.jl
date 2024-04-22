@testset "simple bilevel" begin
    qpn = setup(:simple_bilevel; gen_solution_map=true)

    w1 = [-2.0, -3.0]; x1 = [[-2.0, 0],]
    w2 = [0, -1.0];    x2 = [[0.0, 0],]
    w3 = [1.0, -3.0];  x3 = [[0.0, 0],]
    w4 = [1.0, -1.0];  x4 = [[0.0 ,0],]
    w5 = [1.0, 0.0];   x5 = [[0.5, 0.5],]
    w6 = [0.0, 1.0];   x6 = [[0.5, 0.5], [0.0, 0]]
    w7 = [-1.0, 1+sqrt(2.0)]; x7 = [[-1.0, 0], sqrt(2.0)/2 * ones(2)]
    w8 = [0.0, 0];     x8 = [[0.0, 0],]
    x0 = zeros(2)

    W = [w1, w2, w3, w4, w5, w6, w7, w8]
    X = [x1, x2, x3, x4, x5, x6, x7, x8]
    S = [1, 2, 1, 2, 1, 1, 1, 3]
    for (w, x, s) ∈ zip(W,X,S)
        ret = solve(qpn, [w; x0])
        @test any(isapprox(ret.x_opt, [w; xi]; atol=1e-4) for xi in x)
        @test ret.Sol[2] |> collect |> length ≥ s
    end
end
