using SparseArrays

@testset "convex hull" begin

    ### Test 1 (bounded polyhedra)

    A1 = sparse([1.0 1; 1 -1])
    l1 = fill(1.0, 2)
    u1 = fill(3.0, 2)
    
    A2 = sparse([1.0 0; 0.5 -1; 0.5 1])
    l2 = [-Inf, -2.0, 1]
    u2 = [0.0, Inf, Inf]
    
    A3 = sparse([1.0 0; 0 1])
    l3 = fill(-2.0, 2)
    u3 = fill(-1.0, 2)
    
    P1 = Poly(A1, l1, u1)
    P2 = Poly(A2, l2, u2)
    P3 = Poly(A3, l3, u3)

    PU = QPN.PolyUnion([P1, P2, P3])
    H = QPN.convex_hull(PU)

    (; V, L, R) = QPN.get_verts(H)
    true_verts = [[3.0, 0], [2.0, -1], [-1.0, -2], [-2.0, -2], [-2.0, -1], [-1.0, 1.5], [0.0, 2], [2.0, 1]]
    @test length(V) == 8
    for v in V
        @test any( v ≈ t for t in true_verts)
    end

    ### Test 2 (cones)
    
    A1 = [-1.0 0;
          1 1]
    l1 = [0.0, 0]
    u1 = [Inf, Inf]

    A2 = [1.0 2;
          1.0 -2]
    l2 = l1;
    u2 = u1;

    P1 = Poly(sparse(A1), l1, u1)
    P2 = Poly(sparse(A2), l2, u2)
    PU = QPN.PolyUnion([P1, P2])
    H = QPN.convex_hull(PU) |> QPN.simplify
    (; A, l, u) = QPN.vectorize(H)
    @test any([1.0, 1.0] ≈ r for r in eachrow(A))
    @test any([1.0, 2.0] ≈ r for r in eachrow(A))

    ### Test 3 (rays)
    A1 = [1.0 0] |> sparse
    A2 = [-1.0 0] |> sparse
    l = [0.0,]
    u = [Inf,]
    P1 = Poly(A1, l, u)
    P2 = Poly(A2, l, u)
    PU = QPN.PolyUnion([P1, P2])
    display(PU)
    H = QPN.convex_hull(PU) |> QPN.simplify
    @test isempty(H.poly) # Empty set of slices implies all of ℝ²
end

