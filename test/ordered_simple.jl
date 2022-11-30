@testset "problem of ordered preference simple" begin
    """ 
    QPi: min 0.5*sᵢ^2
        sᵢ, uᵢ
        s.t. -aᵢ'x-bᵢ ≤ aᵢ'(∑ⱼuⱼ) - sᵢ ≤ ∞
        s.t. || ∑ⱼuⱼ ||₁ ≤ 1 (if i == N)


        -inf < a' B*u <= -a'x

    variables (uᵢ sᵢ)ᵢ ∈ 1...N

    """
    N = 4
    allocation = 3.1

    sets = Dict{Int, Poly}()
    qps = Dict{Int, QP}()

    for i in 1:N
        a = [1.0,]
        b = -1.0
        e = zeros(1,N)
        e[i] = 1.0
        A = sparse( [a[1]*e -e] )
        l = [-b,]
        u = [Inf,]
        S = Poly(A, l, u)
        sets[i] = S

        Q = spzeros(2*N, 2*N)
        Q[N+i,N+i] = 1.0
        q = zeros(2*N)
   
        e = zeros(1,N)
        for j = i:N
            e[j] = 1.0
        end
        A = sparse([e zeros(1,N)])
        l = [0.0,]
        u = [allocation,]
        S = Poly(A, l, u)
        sets[N+i] = S
        
        cons = Dict(i=>1.0, N+i=>1.0)
        qps[i] = QP(Quadratic(Q, q), cons, collect(i:N:2N))
    end

    net = [Set([i,]) for i ∈ 1:N]
    qp_net = QPNet(qps, sets, net)

    var_init = zeros(2*N)
    x, Sol = solve(qp_net, var_init; debug=true)
end
