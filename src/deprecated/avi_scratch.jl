
function pivot!(T, row, col)
    pivot = T[row, :] ./ T[row, col]
    T .-= T[:,col]*pivot'
    T[row, :] .= pivot
    return
end

"""
Project x into box (l,u)
"""
function proj(x, l, u)
    @. max(min(x, u), l)
end


function solve_avi(avi::AVI, z0, w; tol=1e-6, max_iters=50)
    q = avi.N*w+avi.o
    n = length(q)
    z = proj(z0, avi.l, avi.u)
    r = avi.M*z + q + z0 - z
    zuvt = zeros(3n+1)
    T = [avi.M -I(n) I(n) r r-q]
    zuvt[1:n] = z0

    basis = []
    basis_lb = zeros(n)
    basis_ub = zeros(n)
    nonbasis = Set([3n+1,])
    for i in 1:n
        if z0[i] - avi.l[i] < tol
            push!(basis, n+i)
            basis_ub[i] = Inf
            push!(nonbasis, i)
            push!(nonbasis, 2n+i)
            zuvt[n+i] = z[i]-z0[i]
        elseif avi.u[i] - z0[i] < tol
            push!(basis, 2n+i)
            basis_ub[i] = Inf
            push!(nonbasis, i)
            push!(nonbasis, n+i)
            zuvt[2n+i] = z0[i] - z[i]
        else
            push!(basis, i)
            basis_lb[i] = avi.l[i]
            basis_ub[i] = avi.u[i]
            push!(nonbasis, n+i)
            push!(nonbasis, 2n+i)
        end
    end
    T .= collect(T[:, basis]) \ T

    entering_ind = 3n+1
    entering_val = 0.0
    entering_lb = -Inf
    entering_ub = 1.0
    entering_dir = 1.0
    
    for iter in 1:max_iters
        delete!(nonbasis, entering_ind)
        d = entering_dir * T[:, entering_ind]
        rhs = T[:, end] - T[:, collect(nonbasis)]*zuvt[collect(nonbasis)] - T[:,entering_ind]*entering_val
        max_val = zeros(n)
        boundary = zeros(n)
        for i in 1:n
            if d[i] > tol
                max_val[i] = (rhs[i]-basis_lb[i]) / d[i]
                boundary[i] = -1
            elseif d[i] < -tol
                max_val[i] = (rhs[i]-basis_ub[i]) / d[i]
                boundary[i] = 1
            else
                max_val[i] = Inf
            end
        end
        min_ind = argmin(max_val)
        min_val = max_val[min_ind]
        if iter == 1 && min_val > Inf 
            zuvt[entering_ind] = 1
            zuvt[basis,1] = T[:, end] - T[:, collect(nonbasis)]*zuvt[collect(nonbasis)]
            return (; z = zuvt[1:n], u=zuvt[n+1:2n], v = zuvt[2n+1:3n], status=SUCCESS)
        elseif isinf(min_val)
            return (; z = zuvt[1:n], u=zuvt[n+1:2n], v = zuvt[2n+1:3n], status=RAY_TERM)
        elseif min_val > (entering_ub-entering_lb)
            exiting_ind = entering_ind
            if entering_dir < 1
                zuvt[exiting_ind] = entering_lb
            else
                zuvt[exiting_ind] = entering_ub
            end
        else
            pivot!(T, min_ind, entering_ind)
            exiting_ind = basis[min_ind]
            if boundary[min_ind] < 0
                zuvt[exiting_ind] = basis_lb[min_ind]
            else
                zuvt[exiting_ind] = basis_ub[min_ind]
            end
            basis[min_ind] = entering_ind
            basis_lb[min_ind] = entering_lb
            basis_ub[min_ind] = entering_ub
        end
        push!(nonbasis, exiting_ind)
        zuvt[basis, 1] = T[:, end] - T[:,collect(nonbasis)] * zuvt[collect(nonbasis)]
        if exiting_ind == 3*n+1
            return (; z = zuvt[1:n], u=zuvt[n+1:2n], v = zuvt[2n+1:3n], status=RAY_TERM)
        elseif exiting_ind > 2*n
            entering_ind = exiting_ind - 2*n
            entering_val = 0.0
            entering_dir = 1.0
            entering_ub = Inf
            entering_lb = 0.0
        elseif exiting_ind > n
            entering_ind = exiting_ind - n
            entering_val = 0.0
            entering_dir = 1.0
            entering_ub = Inf
            entering_lb = 0.0
        elseif boundary[min_ind] < 0
            entering_ind = exiting_ind + n
            entering_val = avi.l[min_ind]
            entering_dir = 1.0
            entering_ub = avi.u[min_ind]
            entering_lb = avi.l[min_ind]
        else
            entering_ind = exiting_ind + 2*n
            entering_val = avi.u[min_ind]
            entering_dir = -1.0
            entering_ub = avi.u[min_ind]
            entering_lb = avi.l[min_ind]
        end
    end 
    return (; z = zuvt[1:n], u=zuvt[n+1:2n], v = zuvt[2n+1:3n], status=MAX_ITERS)
end

