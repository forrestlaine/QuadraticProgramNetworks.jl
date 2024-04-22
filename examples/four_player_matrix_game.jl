"""
Four players. Each player:

min { xi ∈ ℝ² } : ∑ⱼ∑ₖ₌ⱼ₊₁ xj' Aᵢⱼₖ xk
"""
function setup(::Val{:four_player_matrix_game}; edge_list=[], seed=2, show_constellations=false, x_overlay=nothing, kwargs...)

    # Seed 2 gives nice difference between Nash (edge_list=[]) and parallel
    # bilevel (edge_list=[(1,2), (3,4)]) settings.

    rng = MersenneTwister(seed)
    
    x1 = QPNets.variables(:x1, 1:2)
    x2 = QPNets.variables(:x2, 1:2)
    x3 = QPNets.variables(:x3, 1:2)
    x4 = QPNets.variables(:x4, 1:2)
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

    constellations = Dict(i=>Dict(j=>randn(rng, 2) for j in 1:4) for i in 1:4)
    #@infiltrate all(all(norm(constellations[i][j]) ≥ 0.9 for j in setdiff(1:4, i)) for i in 1:4)

    if show_constellations
        f = Figure()
        ax1 = Axis(f[1,1], aspect = 1, xgridvisible=false, ygridvisible=false)
        #ax2 = Axis(f[1,3])
        xlims!(ax1, -2.1, 2.1)
        ylims!(ax1, -1.1, 3.1)
        #xlims!(ax2, -5, 5)
        #ylims!(ax2, -5, 5)
        colors = [:blue, :red, :green, :orange]
        markers = [:utriangle, :rtriangle, :dtriangle, :ltriangle]
        if !isnothing(x_overlay)
            new_colors = [RGBf(rand(3)...) for _ in 1:length(x_overlay)]
        end
        for i = 1:4
            scatter!(ax1, constellations[i][i]..., color=colors[i], marker=markers[i], markersize=15)
            for j in setdiff(1:4, i)
                scatter!(ax1, (constellations[i][j] .+ constellations[i][i])..., color=colors[i], marker=markers[j], markersize=15)
                start = [constellations[i][i]...]
                stop = start .+ [constellations[i][j]...]
                len = norm(stop-start)
                frac = 0.2 / len
                lines!(ax1, 
                       [(1-frac)*start[1] + frac*stop[1], frac*start[1]+(1-frac)*stop[1]],
                       [(1-frac)*start[2] + frac*stop[2], frac*start[2]+(1-frac)*stop[2]],
                       color=colors[i], linewidth=1)
            end

            if !isnothing(x_overlay)
                for (k,(edge_list,x)) in enumerate(x_overlay)
                    scatter!(ax1, x[1], x[2],color=new_colors[k],  marker=markers[1], markersize=15)
                    scatter!(ax1, x[3], x[4],color=new_colors[k],  marker=markers[2], markersize=15)
                    scatter!(ax1, x[5], x[6],color=new_colors[k],  marker=markers[3], markersize=15)
                    scatter!(ax1, x[7], x[8],color=new_colors[k],  marker=markers[4], markersize=15)
                end
            end
        end


        elem_1 = MarkerElement(marker = markers[1], markersize = 15, color=:black)
        elem_2 = MarkerElement(marker = markers[2], markersize = 15, color=:black)
        elem_3 = MarkerElement(marker = markers[3], markersize = 15, color=:black)
        elem_4 = MarkerElement(marker = markers[4], markersize = 15, color=:black)
        elem_5 = LineElement(color = colors[1])
        elem_6 = LineElement(color = colors[2])
        elem_7 = LineElement(color = colors[3])
        elem_8 = LineElement(color = colors[4])
        
        Legend(f[1, 2],
            [
             elem_5, elem_6, elem_7, elem_8,
             elem_1, elem_2, elem_3, elem_4,
            ],
            [
            L"Target constellation for node $1$",
            L"Target constellation for node $2$",
            L"Target constellation for node $3$",
            L"Target constellation for node $4$",
            L"Decision variables indexed by $J^1$",
            L"Decision variables indexed by $J^2$",
            L"Decision variables indexed by $J^3$",
            L"Decision variables indexed by $J^4$",
           ],
            #patchsize = (35, 35), 
            rowgap = 10)

        if !isnothing(x_overlay)
            elems_eq = [MarkerElement(marker=:diamond, markersize=15, color=new_colors[i]) for i in 1:length(x_overlay)]
            labels_eq = []
            for (e, (edge_list, x)) in enumerate(x_overlay)
                if isempty(edge_list)
                    push!(labels_eq, "[]")
                else
                    push!(labels_eq, string(edge_list))
                end
            end
            Legend(f[1,3],
                   elems_eq,
                   labels_eq,
                   patchsize = (35,35), rowgap = 10)
        end

        save("constellations.png", f, px_per_unit = 100)
        display(f)
    end

    for i = 1:4

        #cons = [x[i]; sum(x[i])]
        #lb = [0.0, 0.0, -Inf]
        #ub = [Inf, Inf, 1.0]
        cons = x[i]
        lb = 5*[-1.0, -1.0]
        ub = 5*[1.0, 1.0]
        con_id = QPNets.add_constraint!(qp_net, cons, lb, ub)

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

        QPNets.add_qp!(qp_net, cost, [con_id], dvars)
    end

    QPNets.add_edges!(qp_net, edge_list)
    QPNets.assign_constraint_groups!(qp_net)
    QPNets.set_options!(qp_net; kwargs...)
    qp_net.default_initialization = zeros(8)
    
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

function vis_equilibria(; seed=495)
    edge_list_ps_unique = compute_unique_edge_lists()
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
    x_opt_inds = [22, 23, 29, 40, 1]
    qpn = setup(:four_player_matrix_game; seed, show_constellations=true)
    #qpn = setup(:four_player_matrix_game; seed, show_constellations=true, x_overlay=zip(edge_list_ps_unique[x_opt_inds], x_opts[x_opt_inds]))
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
                    display([I edge_list_ps_unique[I] avg_costs[i][I] 1.96*sqrt.(m2_costs[i][I])./num_success])
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
    return Vector{Tuple{Int64,Int64}}[ [],
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
        qpn = setup(:four_player_matrix_game; edge_list)
       
        locs = Dict()
        D = length(qpn.network_depth_map)
        spacing = 1.0
        pos = []
        preference = nothing
        v1 = qpn.network_depth_map[1]
        weird = false
        weird_node = 0
        if D == 3
            for v1i in v1
                for ci in qpn.network_edges[v1i]
                    if ci ∈ qpn.network_depth_map[3]
                        weird_node = v1i
                        weird = true
                    end
                end
            end
        end

        for k in 1:D
            v = qpn.network_depth_map[k]
            n = length(v)
            if weird
                shift = 0
            else
                shift = -spacing * (n-1)/2.0
            end
            vs = sort(collect(v))
            if !isnothing(preference) && preference ∈ vs
                ind = findfirst(vs .== preference)
                vs = [preference; vs[1:ind-1]; vs[ind+1:end]]
            end
            if k == 1 && weird
                if vs[1] == weird_node
                    vs = [vs[2:end]; vs[1]]
                end
            end

            for (i,vi) in enumerate(vs)
                locs[vi] = [shift + (i-1)*spacing, -(k-1)]
                if D > 2
                    locs[vi] = [-locs[vi][2], locs[vi][1]]
                end
            end
            if length(qpn.network_edges[vs[1]]) > 0
                preference = first(qpn.network_edges[vs[1]])
            else
                preference = nothing
            end
        end
        cmd = "\\expandafter\\newcommand\\csname graph$idx\\endcsname{%\n\\begin{tikzpicture}[scale=0.5,baseline=(current bounding box.center)]\n\\node[circle,draw,fill=yellow!70] (1) at ($(locs[1][1]),$(locs[1][2])) {};\n"
        for i in 2:4
            cmd *= "\\node[circle, draw] ($i) at ($(locs[i][1]),$(locs[i][2])) {};\n"
        end
        cmd *= "\\graph {"
        for e in edge_list
            cmd *= "($(e[1])) -> ($(e[2])); "
        end
        cmd *= "};\n\\end{tikzpicture}\n}\n"
        print(cmd)


        #\newcommand{\graph_22}{%
        #    \begin{tikzpicture}[scale=0.65]
        #        \node[circle, draw, fill=yellow!70] (1) at (0, -1) {};
        #        \node[circle, draw] (2) at (0,-2) {};
        #        \node[circle, draw] (3) at (0,-3) {};
        #        \node[circle, draw] (4) at (0,-4) {};
        #        \graph { (1) -> (2) -> (3) -> (4)};
        #    \end{tikzpicture}
        #}
        

        #for e in edge_list
        #    add_edge!(g, e[1], e[2])
        #end
        #t = TikzGraphs.plot(g, ["α", "β", "γ", "α"], node_style="draw, rounded corners, fill=blue!10", node_styles=Dict(1=>"fill=yellow!10"))
        

        #t = TikzGraphs.plot(g, ["", "", "", ""], graph_options="nodes={draw,circle}", edge_style="line width=1.5pt", node_style="draw", node_styles=Dict(1=>"fill=yellow!20"))
        #t = TikzGraphs.plot(g, ["", "", "", ""], graph_options="nodes={draw,circle}, grow down, branch right", node_style="draw", node_styles=Dict(1=>"fill=yellow!70"))

        #ind = findfirst("layered layout,", t.data)
        #t.data = t.data[1:ind[1]-1] * t.data[ind[end]+1:end]
        
        #same_layer_string = "\n"
        #for (k,v) in qpn.network_depth_map
        #    length(v) < 2 && continue
        #    #k > 1 && continue
        #    same_layer_string *= "{ [same layer] "
        #    for vi in v
        #        same_layer_string *= "$vi, "
        #    end
        #    same_layer_string = same_layer_string[1:end-2]
        #    same_layer_string *= "};\n"
        #end
        #ind = findfirst(";", t.data)[1]
        #data = t.data[1:ind] * same_layer_string * t.data[ind+1:end]
        #t.data = data
        #print(t.data)
        #TikzPictures.save(PDF("graph_$idx"), t)
    end
end
