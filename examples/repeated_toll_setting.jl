"""
A reactive toll setting problem.
N = number of toll lane entrances/exits
L = number of lanes of traffic (L-1 are toll lanes)
0 ≤ p0 ≤ 1 : fraction of traffic in initial toll lane
0 ≤ p[n] ≤ 1 : fraction of traffic in nth toll lane
0 ≤ T[n] : the toll price (in standardized units) for traveling in the nth toll lane
0 ≤ α : the cost factor for lane change congestion between toll sections
0 ≤ D : the extra price (in standardized units) for traveling in any non-toll lane

"""
function setup(::Val{:repeated_toll_setting}; N=2, L = 3, D=0.0, α=0.1, kwargs...)
    T = QPN.variables(:T, 1:L-1, 1:N)
    P = QPN.variables(:P, 1:L, 1:N)
    P0 = ones(L) ./ L

    qp_net = QPNet(T, P)
   
    running_costs = []
    running_revenues = []

    p_players = Dict()
    t_players = Dict()
    for n = N:-1:1
        pnm1 = n > 1 ? P[:, n-1] : P0
        cons = [P[:, n]; sum(P[:,n])]
        lb = [zeros(L); 1.0]
        ub = [fill(Inf, L); 1.0]
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        lane_costs = [T[:,n]; D] + P[:,n]
        lane_change = (P[:,n]-pnm1)'*(P[:,n]-pnm1)
        cost = P[:,n]'*lane_costs + α*lane_change
        push!(running_costs, cost)
        p_players[n] = QPN.add_qp!(qp_net, sum(running_costs), [con_id,], P[:,n])

        revenue = -T[:,n]'*P[1:L-1,n]
        cons = T[:,n]
        lb = fill(0.0, L-1)
        ub = fill(Inf, L-1)
        push!(running_revenues, revenue)
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
        t_players[n] = QPN.add_qp!(qp_net, sum(running_revenues), [con_id,], T[:,n])
    end

    edge_list = []
    for n in 1:N
        push!(edge_list, (t_players[n], p_players[n]))
        if n < N
            push!(edge_list, (p_players[n], t_players[n+1]))
        end
    end
    
    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; exploration_vertices=10)
    qp_net.default_initialization = [zeros((L-1)*N); fill(1.0/L, L*N)]

    qp_net
end
