@testset "my example 2" begin
    """
    variables := [g1 g2 x1 x2 {[zi1 zi2 zi3 zi4 yi]} i ∈ 1...8]
    R1 = { 0 ≤ x1 ≤ 1, 0 ≤ x2 ≤ 1 }
    R2 = { 1 ≤ x1 ≤ 2, 0 ≤ x2 ≤ 1 }
    R3 = { 2 ≤ x1 ≤ 3, 0 ≤ x2 ≤ 1 }
    R4 = { 0 ≤ x1 ≤ 1, 1 ≤ x2 ≤ 2 }
    R5 = { 2 ≤ x1 ≤ 3, 1 ≤ x2 ≤ 2 }
    R6 = { 0 ≤ x1 ≤ 1, 2 ≤ x2 ≤ 3 }
    R7 = { 1 ≤ x1 ≤ 2, 2 ≤ x2 ≤ 3 }
    R8 = { 2 ≤ x1 ≤ 3, 2 ≤ x2 ≤ 3 }
    Ri = { l1i ≤ x1 ≤ u1i, l2i ≤ x2 ≤ u2i }
   
    [R6 R7 R8
     R4 xx R5
     R1 R2 R3]

     fᵢ₁ = min_zi1 zi1*(l1i-x1)
           st      0 ≤ zi1 ≤ 1
     fᵢ₂ = min_zi2 zi2*(x1-u1i)
           st      0 ≤ zi2 ≤ 1
     fᵢ₃ = min_zi3 zi3*(l2i-x2)
           st      0 ≤ zi3 ≤ 1
     fᵢ₄ = min_zi4 zi4*(x2-u2i)
           st      0 ≤ zi4 ≤ 1
     fᵢ₅ = min_yi -yi
           st      y_i <= zij, j∈{1,2,3,4}
     f   = min_{x1,x2}  0.5((x1-g1)^2 + (x2-g2)^2)
           st ∑ᵢ yi >= 1 
    """
    num_vars = 4 + 5*8
    sets = Dict{Int, Poly}()
    qps = Dict{Int, QP}()
    Q = sparse([1 0 -1 0.0;
                0 1 0 -1;
                -1 0 1 0;
                0 -1 0 1])
    Q = [Q spzeros(4, 5*8);
         spzeros(5*8, num_vars)]
    q = zeros(num_vars)
    f = Quadratic(Q, q)
    set_id = 1 
    qp_id = 1
    a = spzeros(1, num_vars)
    for i in 1:8
        a[1,4+i*5] = 1.0
    end
    S = Poly(a, [0.5,], [Inf,])
    sets[set_id] = S
    set_id += 1

    qp = QP(f,Dict(set_id=>1.0),[3,4])
    qps[qp_id] = qp
    qp_id += 1

    l1 = [0.0, 1, 2, 0, 2, 0, 1, 2]
    u1 = l1 .+ 1.0
    l2 = [0.0, 0, 0, 1, 1, 2, 2, 2]
    u2 = l2 .+ 1.0
    b = [l1, u1, l2, u2]
    ϵ = 1e-5
    for i ∈ 1:8
        for (xid, offset, sign, bound) ∈ [(3, 1, -1.0, l1[i]), (3, 2, 1.0, -u1[i]), (4, 3, -1.0, l2[i]), (4, 4, 1.0, -u2[i])]
            ind = 4+(i-1)*5+offset
            Q = spzeros(num_vars, num_vars)
            Q[xid,ind] = sign
            Q[ind,xid] = sign
            q = zeros(num_vars)
            q[ind] = bound-ϵ
            a = spzeros(1,num_vars)
            a[ind] = 1.0
            S = Poly(a, [0.0,], [1.0,])
            sets[set_id] = S
            qp = QP(Quadratic(Q,q), Dict(set_id=>1.0), [ind,])
            qps[qp_id] = qp
            qp_id += 1
            set_id += 1
        end
        ind = 4+i*5
        Q = spzeros(num_vars, num_vars)
        q = zeros(num_vars)
        q[ind] = -1.0
        a = spzeros(4, num_vars)
        for j in 1:4
            a[j, ind-j] = 1.0
            a[j, ind] = -1.0
        end
        u = fill(Inf, 4)
        l = zeros(4)
        S = Poly(a, l, u)
        sets[set_id] = S
        qp = QP(Quadratic(Q,q), Dict(set_id=>1.0), [ind,])
        qps[qp_id] = qp
        qp_id += 1
        set_id += 1
    end
    
    network = [Set([1,]), Set(1+1:1+5*8)]
    qep = QPNet(qps, sets, network)
    x = zeros(num_vars)
    x[1:2] = [2.5,2.6]
    x[3:4] = [0.1, 0.2]
    x, Sol = solve(qep, x; debug=true)
    display(x[3:4])
end
