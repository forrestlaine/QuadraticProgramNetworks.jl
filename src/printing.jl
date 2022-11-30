function display_debug(level, iter, x, num_low, num_tossed)
    space = "    "^level
    if isnothing(num_low) || isnothing(num_tossed)
        println(space*"Level: ", level, ". Iteration: ", iter, ".  x after is: ")
    else
        println(space*"Level: ", level, ". Iteration: ", iter, ". Reasoned about ", num_low, " pieces. Disregarded ", num_tossed, ". x after is: ")
    end
    str = space*"["
    for xi in x
        str *= "%5.2f "
    end
    str *= "]\n"
    format = Printf.Format(str)
    Printf.format(stdout, format, x...)
end

function Base.show(io::IO, ::MIME"text/plain", poly::Poly)
    indent=get(io, :indent, 0)
    space = " "^indent
    n = length(poly)
    d = embedded_dim(poly)
    half = Int(ceil(n/2))
    println(io, space*"Polyhedron in ℝ^", d, " with ", n, " constraints.")
    if n ≤ 50 && d ≤ 50
        (A,l,u,rl,ru) = vectorize(poly)
        for i = 1:n
            str = space*"%5.2f %2s | "
            args = [l[i], rl[i]]
            for j in 1:d
                if iszero(A[i,j])
                    str *= "  ⋅   "
                else
                    str *= "%5.2f "
                    push!(args, A[i,j])
                end
            end
            if i == half
                str *= "| x %2s %5.2f"
            else
                str *= "|   %2s %5.2f"
            end
            push!(args, ru[i])
            push!(args, u[i])
            str *= "\n"
            format = Printf.Format(str)
            Printf.format(io, format, args...)
        end
    end
end

function Base.show(io::IO, mime::MIME"text/plain", node::IntersectionRoot)
    indent=get(io, :indent, 0)
    space = " "^indent
    println(io, space*"Intersection root with ", potential_length(node), " potential polys") 
    for (e,child) in enumerate(node.children)
        ioc = IOContext(io, :indent=>indent+2, :nest=>e==length(node.children))
        show(ioc, mime, child)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", node::IntersectionNode)
    indent=get(io, :indent, 0)
    nest=get(io, :nest, false)
    show(io, mime, node.contributing_poly)
    if nest
        ioc = IOContext(io, :indent=>indent+2)
        for (e,child) in enumerate(node.children)
            ioc = IOContext(io, :indent=>indent+2, :nest=>e==length(node.children))
            show(ioc, mime, child)
        end
    end
end

function Base.show(io::IO, mime::MIME"text/plain", pu::PolyUnion)
    println(io, "Union of ", length(pu), " polyhedra:")
    ioc = IOContext(io, :indent=>2)
    for p in pu
        show(ioc, mime, p)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", f::Quadratic)
    indent=get(io, :indent, 0)
    space = " "^indent
    println(io, space*"Quadratic function:")
    n = length(f.q)
    half = Int(ceil(n/2))

    for i = 1:n
        args = []
        if i == half
            str = space*"0.5 x' |"
        else
            str = space*"       |"
        end
        for j in 1:n
            if iszero(f.Q[i,j])
                str *= "  ⋅   "
            else
                str *= "%5.2f "
                push!(args, f.Q[i,j])
            end
        end
        if i == half
            str *= "| x + x' |"
        else
            str *= "|        |"
        end

        str *= "%5.2f|"
        push!(args, f.q[i])
        str *= "\n"
        format = Printf.Format(str)
        Printf.format(io, format, args...)
    end
end

function Base.show(io::IO, mime::MIME"text/plain", qp::QP)
    indent=get(io, :indent, 0)
    space = " "^indent
    space1 = space*"    "
    println(io, space*"Quadratic program:")
    println(io, space*"__________________")
    println(io, space1*"min  f(x)")
    println(io, space1*"xᵢ : i ∈ ℐ")
    println(io, space1*"s.t. x ∈ S")
    println(io, space*"__________________")
    ioc = IOContext(io, :indent=>indent+2)
    println(io, space*"function f:")
    show(ioc, mime, qp.f)
    println(io, space*"Set S:")
    show(ioc, mime, qp.S)
    println(io,"")
    println(io, space*"Decision variable index set ℐ:")
    print(io, space*"  {")
    for i in qp.indices
        print(io, i, ", ")
    end
    println(io, "}")
end
