#=
The following code contains an implementation of 
Newton-Raphson method to find the root of a 
vector valued function of range and domain n

Algorithm:
    1. Start with an initial guess x0(k=0)
    2. Solve the linear equation AΔx = -f(x(k)), where A 
        is the Jacobian at x0
    3. update x(k+1) = x(k) + Δx
    4. Repeat 2 and 3 until convergence 
=#
using LinearAlgebra
using ForwardDiff

function NewtonRaphson(func, x0; tol=1e-9, MAX_ITERS=50)
    k = 0
    xk = x0 
    
    while k < MAX_ITERS 
        f_xk = func(xk)
        if norm(f_xk) < tol 
            return xk, k 
        end
            
        A = ForwardDiff.jacobian(func, xk)
        if abs(det(A)) < tol 
            println("Newton-Raphson method did not converge; Jacobian is singular")
            return xk, k 
        end 
        Δx = -A \ f_xk # Solving A * Δx = -f_xk 
        
        xk = xk + Δx
        k += 1
    end

    println("Newton-Raphson method did not converge; Maximum Iterations reached") 
    return xk, k 
end

function main()
    # Define the non-linear system of equations
    function f(x) # f:R^4 -> R^4
        f1 = x[1]+2*x[2]-x[1]*(x[1]+4*x[2])-x[2]*(4*x[1]+10*x[2])+3;
        f2 = 3*x[1]+4*x[2]-x[1]*(x[1]+4*x[2])-x[2]*(4*x[1]+10*x[2])+4;
        f3 = 0.5*cos(x[1])+x[3]-sin(x[3])^7;
        f4 = -2*x[2]^2*sin(x[1])+x[4]^3;
        return [f1; f2; f3; f4]
    end

    # Initial guess
    x0 = [-2.; 3.; π; -1.]
    
    root, iterations = NewtonRaphson(f, x0)
    println("Root: $root, iterations: $iterations")
end

# call the main function
main()

