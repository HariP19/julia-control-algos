#=
The following code contains two iterative methods for root finding for a 
    non-linear function f(x), where f(x) is a scalar real-valued function
=#    
using ForwardDiff

# Bisection algo for root finding
function Bisection(f, a, b; θ=1e-9, MAX_ITERS=100)
    if f(a)*f(b) > 0
        error("The function must have different signs at a and b")
    end 

    N = 0 
    
    while(N < MAX_ITERS)
        c = (a + b) / 2
        f_c  = f(c)
        
        # Check for convergence
        if abs(f_c) < θ 
            return c, N
        end
        
        # Check where the root lies betwee and update a or b
        if(sign(f_c)* sign(f(a)) > 0)
            a = c 
        else
            b = c 
        end
       
        N += 1
    end

    println("Max iterations reached, solution not found")
    return c, N
end

# Newton algo for root finding
function Newton(f, x_0; θ = 1e-9, MAX_ITERS=100)
    x_i = x_0
    N = 0
    
    while N < MAX_ITERS
        f_i = f(x_i)
        df_i = ForwardDiff.derivative(f,x_i)
        
        if df_i == 0
            println("derivative is zero. No solution found.")
            return x_i, N
        
        x_n = x_i - f_i / df_i
        
        # Check for convergence
        if abs(f_i)<θ
            return x_i, N
        end

        # Update x_i
        x_i = x_n 
        N += 1
    end

    println("Max iterations reached, solution not found")
    return x_i, N
end

function  main()
    # non-linear function you want to solve
    function f(x)
        return (x-1)^2 - 4
    end 

    # Bisection solution
    a = 1; b = 30
    root, iterations = Bisection(f, a, b)
    println("Bisection: Root: $root, Iterations: $iterations")

    # Newton solution 
    x_0 = 7
    root, iterations = Newton(f, x_0)
    println("Newton: Root: $root, Iterations: $iterations")
    
end

# Call the main function
main()
