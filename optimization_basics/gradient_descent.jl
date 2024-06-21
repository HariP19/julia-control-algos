#=
The following code contains algorithmic implementation of 
Gradient Descent for minimizing an objective function 
    f: R^m -> R 

Algorithm:
    1. Start with an initial guess x0 (k = 0)
    2. Evaluate ∇f(xk). If ||∇f(xk)|| = 0, then the algorithm is converged 
    3. Update the decision variable via x(k+1) = xk - α∇f(xk)
    4. Repeat until convergence  
=#
using ForwardDiff
using LinearAlgebra
using Plots 

function GradientDescent(func, x0; α=0.1, MAX_ITERS=100, tol=1e-9)
    xk = x0 
    k = 0

    while k < MAX_ITERS
        ∇f_k = ForwardDiff.gradient(func, xk) 
        
        # Check convergence using the norm of gradient
        if norm(∇f_k) < tol
            display(∇f_k)
            return xk, k 
        end
        
        # Update the parameter vector
        xk = xk - α * ∇f_k
        k += 1
    end
    println("Convergence failed: Max iterations reached")
    return xk, k
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

    x0 = [1; 0]
    x_arg_min, iterations = GradientDescent(f2, x0, α=0.001, MAX_ITERS=1000)

    println("Arg min of f1(x): $x_arg_min")
    println("Number of iterations: $iterations")

end

# Call main function
main()


