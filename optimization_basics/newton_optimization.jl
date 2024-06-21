#= 
This code contains the algorithm for Newton optimization, which 
    finds the parameter vector which corresponds to the min / max
    of an objective function. 

Algorithm: 
        1. Start with an initial guess x₀ (k = 0)
        2. If ||∇f(xₖ)|| = 0, then the algorithm is converged. 
        3. Solve the linear equation HₖΔx = -∇f(xₖ).
        4. Update xₖ to xₖ + α*Δx 
        5. Repeat until convergence. 
=# 
using ForwardDiff
using LinearAlgebra

function NewtonOptimization(func, x0; α=0.99, tol=1e-9, MAX_ITERS=100)
    # Initialization
    xk = x0 
    k = 0

    while k < MAX_ITERS
        # Evaluate gradient and hessian at xk 
        ∇f_k = ForwardDiff.gradient(func, xk)
        H_k = ForwardDiff.hessian(func, xk)
        
        # Check for convergence 
        if norm(∇f_k) < tol
            if abs(det(H_k)) > tol 
                if all(sign.(eigen(H_k).values) .> 0)
                    println("Converged to a local minimum at iteration: k = $k")
                    return xk
                elseif all(sign.(eigen(H_k).values) .< 0) 
                    println("Converged to a local maximum at iteration: k = $k")
                    return xk
                else
                    println("Converged to a saddle point at iteration: k = $k")
                    return xk
                end 
            else 
                println(" Insufficient information to proceed. Newton-Raphson method did not converge; Hessian is singular.")
                return xk
            end
        end 

        # Update the parameter vector
        dx = -H_k \  ∇f_k  #Solve HₖΔx = -∇f(xₖ)
        xk = xk + α * dx 
        k += 1
    end 
    println("Convergence failed: Max iterations reached")
    return xk
end

function main()
    # Simple objective function 
    function f1(x)
        @assert length(x) == 2
        A = [2 -1; -1 2]
        return dot(x, A * x)
    end

    # Complex objective function
    function f2(x)
        @assert length(x) == 2
        x1, x2 = x
        return -(x1^4 - x2^4) * exp(-0.1*x1^2 - 0.1*x2^2)
    end

    x0 = [3.5; 0]
    x_arg_min = NewtonOptimization(f2, x0)

    println("Arg min of f1(x): $x_arg_min")

end

main()