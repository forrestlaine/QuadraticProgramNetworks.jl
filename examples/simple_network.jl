"""
variables := [x1, x2, x3]
f1: (x1)² + (x2-1)²
f2: (x2+1)²
f3: (x3)²
C3: { x : x2 - x1 - x3 ≥ 0 }
"""
function setup(::Val{:simple_network}; edge_version=1, kwargs...)
    
    x = Symbolics.variables(:x, 1:3)
    
    qp_net = QPNet(x)
  
    lb = [0.0,]
    ub = [Inf,]
    cons = [x[2]-x[1]-x[3],]
    con_id = QPN.add_constraint!(qp_net, cons, lb, ub)
           
    cost_1 = (x[1])^2 + (x[2]-1)^2
    dvars_1 = [x[1],]
    con_ids_1 = []
    qp_id1 = QPN.add_qp!(qp_net, cost_1, con_ids_1, dvars_1...)
    
    cost_2 = (x[2]+1)^2
    dvars_2 = [x[2],]
    con_ids_2 = []
    qp_id2 = QPN.add_qp!(qp_net, cost_2, con_ids_2, dvars_2...)
    
    cost_3 = (x[3])^2
    dvars_3 = [x[2],x[3]]
    con_ids_3 = [con_id,]
    qp_id3 = QPN.add_qp!(qp_net, cost_3, con_ids_3, dvars_3...)

    edge_list_v1 = [(qp_id2,qp_id3)]
    edge_list_v2 = [(qp_id1,qp_id3),(qp_id2,qp_id3)]
    edge_list_v3 = [(qp_id1,qp_id2),(qp_id2,qp_id3)]
    edge_versions = Dict(1=>edge_list_v1, 2=>edge_list_v2, 3=>edge_list_v3)

    QPN.add_edges!(qp_net, edge_versions[edge_version])
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end

