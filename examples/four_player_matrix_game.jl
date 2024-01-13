"""
Four players. Each player:

min { xi ∈ ℝ² } : ∑ⱼ∑ₖ₌ⱼ₊₁ xj' Aᵢⱼₖ xk
"""
function setup(::Val{:four_player_matrix_game}; edge_list=[], seed=2, show_constellations=false, kwargs...)

    # Seed 2 gives nice difference between Nash (edge_list=[]) and parallel
    # bilevel (edge_list=[(1,2), (3,4)]) settings.

    rng = MersenneTwister(seed)
    
    x1 = Symbolics.variables(:x1, 1:2)
    x2 = Symbolics.variables(:x2, 1:2)
    x3 = Symbolics.variables(:x3, 1:2)
    x4 = Symbolics.variables(:x4, 1:2)
    x = Dict(1=>x1, 2=>x2, 3=>x3, 4=>x4)
    
    qp_net = QPNet(x1,x2,x3,x4)
  
    A = Dict() 
    for i = 1:4
        for j = 1:4
            for k = j+1:4
                A[i,j,k] = rand(rng,1:4,2,2)
            end
        end
    end

    #locations = [[-0.5,-0.5], [-0.5, 0.5], [0.5, 0.5], [0.5,-0.5]]
    #locations = [randn(rng, 2) for i in 1:4]
    #constellations = Dict(i=>Dict(j=>locations[mod1(j+i-1,4)] for j in 1:4) for i in 1:4)
    constellations = Dict(i=>Dict(j=>randn(rng, 2) for j in 1:4) for i in 1:4)

    #for i = 1:4
    #    constellations[i][i] = zeros(2)
    #end

    if show_constellations
        f = Figure()
        ax1 = Axis(f[1,1])
        ax2 = Axis(f[1,2])
        ax3 = Axis(f[2,1])
        ax4 = Axis(f[2,2])
        ax = [ax1, ax2, ax3, ax4]
        for i = 1:4
            scatter!(ax[i], constellations[i][1]..., color=:blue)
            scatter!(ax[i], constellations[i][2]..., color=:red)
            scatter!(ax[i], constellations[i][3]..., color=:green)
            scatter!(ax[i], constellations[i][4]..., color=:yellow)
            xlims!(ax[i], -1, 1)
            ylims!(ax[i], -1, 1)
        end
        display(f)
    end

    for i = 1:4

        #cons = [x[i]; sum(x[i])]
        #lb = [0.0, 0.0, -Inf]
        #ub = [Inf, Inf, 1.0]
        cons = x[i]
        lb = 5*[-1.0, -1.0]
        ub = 5*[1.0, 1.0]
        con_id = QPN.add_constraint!(qp_net, cons, lb, ub)

        ent = [0.5, 0.5]
        #cost = 0
        #for j = 1:4
        #    for k = j+1:4
        #        cost += x[j]'*A[i,j,k]*x[k]
        #    end
        #end
        cost = 0
        #targs = map(1:4) do j
        #    a = 0.5*randn(rng, 2)
        #    #a ./ sum(a)
        #    a
        #end

        #for j = 1:4
        #    for k = 1:4
        #        mult = (j==k) ? 1.0 : 0.25
        #        cost += mult*(x[j]-targs[j])'*(x[k]-targs[k])
        #    end
        #end

        cost = 0
        for j = 1:4
            if j == i
                d = x[i] - constellations[i][j]
            else
                d = x[j] - x[i] - constellations[i][j]
            end
            cost += d'*d
        end

        #all_x = vcat(x[1], x[2], x[3], x[4])
        #Q = randn(rng, 8, 8)
        #Q = Q'*Q
        #q = randn(rng, 8)
        #cost = all_x'*Q*all_x + all_x'*q

        dvars = x[i]

        QPN.add_qp!(qp_net, cost, [con_id], dvars)
    end

    QPN.add_edges!(qp_net, edge_list)
    QPN.assign_constraint_groups!(qp_net)
    QPN.set_options!(qp_net; kwargs...)
    
    qp_net
end


function search_for_game(seed_range)
    # 4080 results in 8 for all_edges=[(1,2),(2,3),(3,4)]
    all_edges = [(1,2),(1,3),(1,4),(2,3),(2,4),(3,4)]
    all_edges = [(4,2),(4,3),(4,1),(2,3),(2,1),(3,1)]
    edge_list_ps = powerset(all_edges) |> collect

    el_dict = Dict()
    for (e,edge_list) in enumerate(edge_list_ps) 
        qpn = setup(:four_player_matrix_game; edge_list)
        if qpn.network_edges in keys(el_dict)
            push!(el_dict[qpn.network_edges], e)
        else
            el_dict[qpn.network_edges] = [e]
        end
    end
    unique_edges = []
    for e in 1:length(edge_list_ps)
        unique=true
        for (k,v) in el_dict
            if e in v && length(v) > 1 && e != minimum(v)
                unique=false
            end
        end
        if unique
            push!(unique_edges, e)
        end
    end
    edge_list_ps_unique = edge_list_ps[unique_edges]

    num_unique_equilibria = map(seed_range) do seed
        x_opts = map(edge_list_ps_unique) do edge_list
            qpn = setup(:four_player_matrix_game; edge_list, seed)
            try
                ret = solve(qpn, fill(0.0, 8))
                ret.x_opt
            catch e
                nothing
            end
        end
        failed = isnothing.(x_opts)
        if sum(failed) > 0
            @info "$seed => $([Inf for i in 1:length(x_opts)])"
            0
        else
            equilibria = Dict(i=>[i] for i in 1:length(x_opts))
            for i in 1:length(edge_list_ps_unique)
                qpn = setup(:four_player_matrix_game; edge_list=edge_list_ps_unique[i], seed)
                for j in 1:length(x_opts)
                    j == i && continue
                    try 
                        ret = solve(qpn, x_opts[j])
                        if ret.x_opt ≈ x_opts[j]
                            @debug "$j is an equilibrium for problem $i"
                            push!(equilibria[i], j)
                        end
                    catch e
                        continue
                    end
                end
            end
            if all(length(equilibria[i]) == 1 for i in 1:length(x_opts))
                @info "FOUND ONE! seed=$seed"
            end
            @info "$seed => $(collect(length(equilibria[i]) for i in 1:length(x_opts)))"
            sum(length(equilibria[i]) == 1 for i in 1:length(x_opts))
        end
    end
    ii = argmax(num_unique_equilibria)
    @info "The best seed was $(seed_range[ii]) with $(num_unique_equilibria[ii]) unique equilibria"
end

function analyze_equilibria(seed_range)
    edge_list_ps_unique = compute_unique_edge_lists()
    @infiltrate

    num_success = 0

    avg_costs = map(1:4) do i
        map(edge_list_ps_unique) do el
            0.0
        end
    end
    m2_costs = map(1:4) do i
        map(edge_list_ps_unique) do el
            0.0
        end
    end

    @showprogress for seed in seed_range
        try
            x_opts = map(edge_list_ps_unique) do edge_list
                qpn = setup(:four_player_matrix_game; edge_list, seed)
                try
                    ret = solve(qpn, fill(0.0, 8))
                    ret.x_opt
                catch e
                    if e isa InterruptException
                        return
                    end
                    nothing
                end
            end

            @infiltrate any(isnothing.(x_opts))

            #f = Figure()
            #ax = Axis(f[1, 1])

            #x_show = x_opts[1:1]
            #scatter!(ax, [x[1] for x in x_show], [x[2] for x in x_show], color=:blue)
            #scatter!(ax, [x[3] for x in x_show], [x[4] for x in x_show], color=:red)
            #scatter!(ax, [x[5] for x in x_show], [x[6] for x in x_show], color=:green)
            #scatter!(ax, [x[7] for x in x_show], [x[8] for x in x_show], color=:yellow)
            #xlims!(ax, -1, 1)
            #ylims!(ax, -1, 1)
            #window = display(GLMakie.Screen(), f)
            #display(window)
            
            if any(isnothing.(x_opts))
                error("At least one edge list failed.")
            end
            num_success +=  1

            qpn = setup(:four_player_matrix_game; seed)
            x_nash = x_opts[1]

            for (e, (x, edge_list)) in enumerate(zip(x_opts, edge_list_ps_unique))
                for i in 1:4
                    if isempty(edge_list)
                        f = qpn.qps[i].f(x)
                    else
                        f = qpn.qps[i].f(x) - qpn.qps[i].f(x_nash)
                    end
                    delta = f - avg_costs[i][e]
                    avg_costs[i][e] += delta / num_success
                    delta2 = f - avg_costs[i][e]
                    m2_costs[i][e] += delta * delta2
                end
            end

            if mod(num_success, 10) == 0
                
                @info "***************************"
                @info "STATUS UPDATE: (seed=$seed)"
                @info "***************************"

                for i in 1:1
                    I = sortperm(avg_costs[i])
                    display([edge_list_ps_unique[I] avg_costs[i][I] 1.96*sqrt.(m2_costs[i][I])./num_success])
                end
            end
            if mod(num_success, 1000) == 0
                output = ""
                for (idx, edge_list) in enumerate(edge_list_ps_unique)
                #for (edge_list, costs) in results
                    output *= "\\{"
                    for e in edge_list
                        output *= "($(e[1]), $(e[2])), "
                    end
                    if length(edge_list) > 0
                        output = output[1:end-2]
                    end
                    cs = [round(avg_costs[i][idx]; digits=4) for i in 1:4]
                    cstd = [round(sqrt(m2_costs[i][idx]) / num_success; digits=4) for i in 1:4]
                    output *= "\\} &"
                    for i in 1:4
                        output *= " $(cs[i])±$(cstd[i]) &"
                    end
                    output = output[1:end-1]
                    output *= " \\\\ \n"
                    #output *= "\\} & $(cs[1])±$(cstd[1]) & $(cs[2])±$(cstd[2]) & $(cs[3]) & $(cs[4]) \\\\ \n"
                end
                print(output)
            end


        catch err
            if err isa InterruptException
                return
            end
            @info "Bad seed: $seed"
            continue
        end
    end
    
    @info "Percent successful seeds: $(100*num_success/length(seed_range))"
end

function graph_is_redundant(edge_list, existing_edge_lists)
    perms = [Dict(1=>1, 2=>3, 3=>4, 4=>2), 
             Dict(1=>1, 2=>2, 3=>4, 4=>3),
             Dict(1=>1, 2=>3, 3=>2, 4=>4),
             Dict(1=>1, 2=>4, 3=>3, 4=>2),
             Dict(1=>1, 2=>4, 3=>2, 4=>3),
             Dict(1=>1, 2=>2, 3=>3, 4=>4)]
    for perm in perms
        el = Set((perm[e1], perm[e2]) for (e1, e2) in edge_list)
        if el ∈ existing_edge_lists
            return true
        end
    end
    return false
end

function compute_unique_edge_lists()
    return [ [],
        [(1, 2)],
        [(2, 3)],
        [(2, 1)],
        [(1, 2), (1, 3)],
        [(1, 2), (2, 3)],
        [(1, 2), (3, 1)],
        [(3, 2), (1, 2)],
        [(1, 2), (3, 4)],
        [(2, 4), (2, 3)],
        [(2, 1), (2, 3)],
        [(3, 1), (2, 3)],
        [(3, 4), (2, 3)],
        [(4, 1), (2, 3)],
        [(4, 3), (2, 3)],
        [(3, 1), (2, 1)],
        [(1, 2), (1, 3), (1, 4)],
        [(2, 4), (1, 2), (1, 3)],
        [(1, 2), (4, 1), (1, 3)],
        [(1, 2), (4, 2), (1, 3)],
        [(2, 4), (1, 2), (2, 3)],
        [(1, 2), (3, 4), (2, 3)],
        [(1, 2), (4, 1), (2, 3)],
        [(1, 2), (4, 2), (2, 3)],
        [(1, 2), (4, 3), (2, 3)],
        [(1, 2), (3, 1), (3, 4)],
        [(1, 2), (3, 1), (4, 1)],
        [(1, 2), (3, 1), (4, 2)],
        [(1, 2), (3, 1), (4, 3)],
        [(3, 2), (1, 2), (3, 4)],
        [(3, 2), (1, 2), (4, 2)],
        [(3, 2), (1, 2), (4, 3)],
        [(2, 4), (2, 1), (2, 3)],
        [(2, 4), (3, 1), (2, 3)],
        [(2, 1), (3, 4), (2, 3)],
        [(4, 1), (2, 1), (2, 3)],
        [(4, 2), (2, 1), (2, 3)],
        [(4, 3), (2, 1), (2, 3)],
        [(3, 1), (4, 1), (2, 3)],
        [(3, 1), (4, 2), (2, 3)],
        [(3, 1), (4, 3), (2, 3)],
        [(3, 1), (4, 1), (2, 1)],
        [(2, 4), (1, 2), (1, 3), (3, 4)],
        [(1, 2), (4, 2), (1, 3), (4, 3)],
        [(1, 2), (3, 1), (4, 2), (3, 4)],
        [(2, 4), (3, 1), (4, 1), (2, 3)],
        [(4, 1), (4, 3), (2, 1), (2, 3)]]

    all_edges = [(1,2),(1,3),(1,4),(2,3),(2,4),(2,1),(3,1),(3,2),(3,4),(4,1),(4,2),(4,3)]
    edge_list_ps = powerset(all_edges) |> collect
    el_set = Set()
    edge_list_unique = []
    for (e, edge_list) in enumerate(edge_list_ps)
        local qpn, network_edges, edge_list_small
        try
            qpn = setup(:four_player_matrix_game; edge_list)
            network_edges = qpn.network_edges 
            edge_list_small = Set{Tuple{Int64,Int64}}()
            depth = length(qpn.network_depth_map)
            for d in 1:depth
                for (k, v) in network_edges
                    if k in qpn.network_depth_map[d]
                        for vi in v
                            push!(edge_list_small, (k,vi))
                        end
                    end
                end
            end
        catch err
            continue
        end
        if graph_is_redundant(edge_list_small, edge_list_unique)
            continue
        else
            push!(edge_list_unique, edge_list_small)
        end
    end
    map(edge_list_unique) do el
        collect(el)
    end
end


function generate_graph_images()
    edge_list_ps_unique = compute_unique_edge_lists()

    for (idx,edge_list) in enumerate(edge_list_ps_unique)
        g = SimpleDiGraph(UInt8(4))
        for e in edge_list
            add_edge!(g, e[1], e[2])
        end
        #t = TikzGraphs.plot(g, ["α", "β", "γ", "α"], node_style="draw, rounded corners, fill=blue!10", node_styles=Dict(1=>"fill=yellow!10"))

        #t = TikzGraphs.plot(g, ["", "", "", ""], graph_options="nodes={draw,circle}", edge_style="line width=1.5pt", node_style="draw", node_styles=Dict(1=>"fill=yellow!20"))
        t = TikzGraphs.plot(g, ["", "", "", ""], graph_options="nodes={draw,circle}", node_style="draw", node_styles=Dict(1=>"fill=yellow!50"))
        TikzPictures.save(PDF("graph_$idx"), t)
    end
end
