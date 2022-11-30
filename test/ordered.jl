@testset "problem of ordered preference toy" begin
    """ 
    QPi: min 0.5*sᵢ^2
        sᵢ, uᵢ
        s.t. -a'x-1 ≤ aᵢ'(∑ⱼuⱼ) - sᵢ ≤ ∞
        s.t. || ∑ⱼuⱼ ||₁ ≤ 1 (if i == N)


        -inf < a' B*u <= -a'x

    variables (uᵢ sᵢ)ᵢ ∈ 1...N

   _______
   x->    \
           \
            \
            |
            |
            |
    
    """

    N = 3 
    x = [0, 0.0]
    v = [1.5, 0.7]
    #xx = [1.1, 1.4]

    xx = x + v
    
    A = sparse([ones(1,N) ones(1,N) zeros(1,N);
                ones(1,N) -ones(1,N) zeros(1,N)])
    l = [-1, -1.0]
    u = [1, 1.0]
    S0 = Poly(A, l, u)

    sets = Dict(0=>S0)
    qps = Dict{Int, QP}()

    for i in 1:N
        θ = (i-1)*π/(2*(N-1))
        a = -[cos(θ), sin(θ)]
        e = zeros(1,N)
        e[i] = -1.0
        A = sparse( [a[1]*ones(1,N) a[2]*ones(1,N) e] )
        l = [-a'*xx-1.0,]
        u = [Inf,]
        S = Poly(A, l, u)
        sets[i] = S

        # TODO is this set needed? why cycling?
        A = sparse( [e zeros(1,N) zeros(1,N);
                     zeros(1,N) e zeros(1,N)] )
        l = [-1, -1.] 
        u = [1, 1.] 
        S = Poly(A, l, u)
        #sets[N+i] = S

        Q = spzeros(3*N, 3*N)
        Q[2*N+i,2*N+i] = 1.0
        q = zeros(3*N)
        
        #cons = Dict(i=>1.0, N+i=>1.0)
        cons = Dict(i=>1.0)
        if i == N 
            cons[0] = 1.0
        end
        qps[i] = QP(Quadratic(Q, q), cons, collect(i:N:3N))
    end

    net = [Set([i,]) for i ∈ 1:N]
    qp_net = QPNet(qps, sets, net)

    
    var_init = zeros(3*N)
    x, Sol = solve(qp_net, var_init; debug=true)
    #Sol = QPN.PolyUnion(collect(Sol))
end
