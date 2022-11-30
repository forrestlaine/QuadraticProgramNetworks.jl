@testset "pessimistic commitment trilevel" begin
    """
    variables := p, q₁, q₂ ∈ Δ¹ ⊂ ℝ²
    f1: p'Aq₁  (p)
    f2: -p'Aq₁ (q₁)
    f3: p'Bq₂ (q₂)
    """
    A = -[1.0 3; 0 2]
    B = -[1.0 0; 0 1]
    Q1 = sparse([spzeros(2,2) A spzeros(2,2);
                A' spzeros(2,4)
                spzeros(2,6)])
    Q2 = -Q1
    Q3 = sparse([spzeros(2,4) B';
                 spzeros(2,6);
                 B spzeros(2,4)])

    q1 = q2 = q3 = zeros(6)
    f1 = Quadratic(Q1, q1)
    f2 = Quadratic(Q2, q2)
    f3 = Quadratic(Q3, q3)

    A1 = sparse([1 1 0 0 0 0;
                 1 0 0 0 0 0;
                 0 1 0 0 0 0.])
    A2 = sparse([0 0 1 0 -1 0;
                 0 0 0 1 0 -1.])
    A3 = sparse([0 0 0 0 1 1;
                 0 0 0 0 1 0;
                 0 0 0 0 0 1.0])
    l = [1.0,0,0]
    u = [1.0,Inf,Inf]
    S1 = Poly(A1, l, u)
    S3 = Poly(A3, l, u)
    l = u = [0.0,0]
    S2 = Poly(A2, l, u)
    sets = Dict(1=>S1, 2=>S2, 3=>S3)
    qp1 = QP(f1, Dict(1=>1.0), [1,2])
    qp2 = QP(f2, Dict(2=>1.0), [3,4])
    qp3 = QP(f3, Dict(3=>1.0), [5,6])
    qps = Dict(1=>qp1, 2=>qp2, 3=>qp3)
    net = [Set(1), Set(2), Set(3)] 
    qp_net = QPNet(qps, sets, net)
    x = [0, 1.0, 0, 1, 0, 1]
    x, Sol = solve(qp_net, x; debug=true, gen_Sol=true)
    xx = [0.1, 0.9, 0, 1, 0, 1]
    Sol = collect(Sol)
    display(Sol)
end
