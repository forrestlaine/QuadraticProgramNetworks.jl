@testset "my example 2" begin
    """
    variables := [g1 g2 x1 x2 {yi ∈ 1...8} a0...3 b0...3]
    R1 = { 0 ≤ x1 ≤ 1, 0 ≤ x2 ≤ 1 }
    R2 = { 1 ≤ x1 ≤ 2, 0 ≤ x2 ≤ 1 }
 
     R1 R2

     fᵢ = min_aᵢ aᵢ*(i-x1)
           st      0 ≤ aᵢ ≤ 1
     gᵢ = min_bᵢ bᵢ*(i-x2)
           st      0 ≤ bᵢ ≤ 1
     h1 = min_y1  -y1
           st      y1 ≤ a0-a1
                   y1 ≤ b0-b1
     h2 = min_y2  -y2
            st     y2 ≤ a1-a2
                   y2 ≤ b0-b1
     ...


     f   = min_{x1,x2}  0.5((x1-g1)^2 + (x2-g2)^2)
           st ∑ᵢ yi >= 1 
    """
    num_vars = 4 + 2 + 2 + 3
    sets = Dict{Int, Poly}()
    qps = Dict{Int, QP}()
    Q = sparse([1 0 -1 0.0;
                0 1 0 -1;
                -1 0 1 0;
                0 -1 0 1])
    Q = [Q spzeros(4, 7);
         spzeros(7, num_vars)]
    q = zeros(num_vars)
    f = Quadratic(Q, q)
    set_id = 1 
    qp_id = 1
    a = spzeros(1, num_vars)
    a[1,5:6] .= 1.0
    S = Poly(a, [0.25,], [Inf,])
    sets[set_id] = S
    set_id += 1

    qp = QP(f,Dict(set_id=>1.0),[3,4])
    qps[qp_id] = qp
    qp_id += 1

    k = 0
    for j ∈ 1:1
        for  i ∈ 1:2
            k += 1
            ind = 4 + k
            Q = spzeros(num_vars, num_vars)
            q = zeros(num_vars)
            q[ind] = -1.0
            a = spzeros(2, num_vars)
            a[1,ind] = 1.0
            a[1,6+i] = -1.0
            a[1,7+i] = 1.0
            a[2,ind] = 1.0
            a[2,9+j] = -1.0
            a[2,10+j] = 1.0
            S = Poly(a, [-Inf, -Inf], [0, 0.0])
            sets[set_id] = S
            qp = QP(Quadratic(Q,q), Dict(set_id=>1.0), [ind,])
            qps[qp_id] = qp
            qp_id += 1
            set_id += 1
        end
    end

    for i ∈ 0:2
        for (xind, ind_offset) ∈ ( (3, 7), (4, 10) )
            i == 2 && xind == 4 && continue
            ind = ind_offset+i
            Q = spzeros(num_vars, num_vars)
            Q[xind,ind] = -1.0
            Q[ind,xind] = -1.0
            q = zeros(num_vars)
            q[ind] = 1.0*i
            a = spzeros(1, num_vars)
            a[ind] = 1.0
            S = Poly(a, [0.0,], [1.0,])
            sets[set_id] = S
            qp = QP(Quadratic(Q,q), Dict(set_id=>1.0), [ind,])
            qps[qp_id] = qp
            qp_id += 1
            set_id += 1
        end
    end
    
    network = [Set([1,]), Set(2:length(qps))]
    qep = QPNet(qps, sets, network)
    x = randn(num_vars)
    x[1:2] = [1.5,0.4]
    x[3:4] = [0.1, 0.2]
    x, Sol = solve(qep, x; debug=true, high_dim=true)
    display(x[3:4])
end
