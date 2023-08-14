using QPN
using Infiltrator

"""
A reactive toll setting problem.
N = number of toll lane entrances/exits
0 ≤ p0 ≤ 1 : fraction of traffic in initial toll lane
0 ≤ p[n] ≤ 1 : fraction of traffic in nth toll lane
0 ≤ T[n] : the toll price (in standardized units) for traveling in the nth toll lane
0 ≤ α : the cost factor for lane change congestion between toll sections
0 ≤ D : the extra price (in standardized units) for traveling in any non-toll lane

"""
function setup(::Val{:repeated_toll_setting}; N=2, D=0.0, p0=0.5, α=0.1, kwargs...)
    T = QPN.variables(:T, 1:N)
    P = QPN.variables(:P, 1:N)

    qp_net = QPNet(T, P)
   
    running_costs = Any[]
    running_revenues = Any[]
    for n = N:-1:1
        pnm1 = n > 1 ? P[n-1] : p0
        cons = [P[n],]
        lb = [0.0,]
        ub = [1.0,]
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        lane_costs = [T[n]+P[n], D+(1.0-P[n])]
        traffic = [P[n], 1.0-P[n]]
        lane_change = (P[n]-pnm1)^2
        cost = traffic'*lane_costs + α*lane_change
        push!(running_costs, cost)
        level = 2*n
        QPN.add_qp!(qp_net, level, sum(running_costs), [con_id,], P[n])

        revenue = -T[n] * P[n]
        cons = [T[n],]
        lb = [0.0,]
        ub = [Inf,]
        push!(running_revenues, revenue)
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        level = 2*n-1
        QPN.add_qp!(qp_net, level, sum(running_revenues), [con_id,], T[n])
    end

    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end
