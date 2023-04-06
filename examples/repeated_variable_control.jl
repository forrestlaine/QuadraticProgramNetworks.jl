using QPN
using SparseArrays

"""
variables := x ∈ ℝⁿ, s ∈ ℝ²


min_x f(x)
s.t. l ≤ Ax ≤ u

Reformulated as 

min_{x,s} f(x)
s.t. x,s ∈ argmin s's
            s.t. 0 ≤ (Ax-l)ᵢ - s₁ ∀ i 
                 0 ≤ (u-Ax)ᵢ - s₂ ∀ i

"""
function setup(; rng=MersenneTwister(1), n=3, m=2, kwargs...)
 
    x = QPN.variables(:x, 1:n) 
    s = QPN.variables(:s, 1:2)

    Q = randn(rng, n, n)
    Q = Q'*Q
    q = randn(rng, n)

    A = randn(rng, m, n)
    l = fill(-1.0, m)
    u = fill( 1.0, m)
 
    qp_net = QPNet(x,s)
  
    lb = fill(0.0, 2*m)
    ub = fill(Inf, 2*m)
    cons = [[A[i,:]'*x-l[i]-s[1] for i=1:m]; [u[i]-A[i,:]'*x-s[2] for i = 1:m]]
    con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
           
    cost = s'*s
    level = 2
    QPN.add_qp!(qp_net, level, cost, [con_id,], x, s)

    cost = 0.5*x'*Q*x + x'*q
    level = 1
    QPN.add_qp!(qp_net, level, cost, [])
    
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)

    (; qp_net, A, l, u, Q, q)
end

