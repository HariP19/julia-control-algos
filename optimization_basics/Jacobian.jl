using LinearAlgebra
using ForwardDiff

function Jacobian(func, x0, h=0.001)
    #Numerical Jacobian of *f: R^m -> R^n*
    m = length(x0)
    f0 = func(x0)
    n = length(f0)

    if m == 1 # f: R -> R^n
        return (func(x0 + h) - func(x0 - h)) / 2*h 
    else
        Im = Matrix(1.0I, m, m)
        A = zeros(n,m)

        for i = 1:m
            ei = Im[:,i:i]
            A[:,i] = (func(x0 + h*ei) - func(x0 - h*ei)) / (2*h);
        end 
        return A
    end
end

function main()
    function f(x)
        x1, x2, x3 = x 
        return [x1*x2*x3; log(2 + cos(x1)) + x2^x1; x1*x3/(1+ x2^2)]
    end

    x0 = [pi; 1.0; 2.0]
    # Measure time for Jacobian
    println("Timing for Jacobian:")
    @time J_manual = Jacobian(f, x0)

    # Measure time for ForwardDiff.jacobian
    println("Timing for ForwardDiff.jacobian:")
    @time J_forwarddiff = ForwardDiff.jacobian(f, x0)

    # display(J_manual)
    # display(J_forwarddiff)
end

main()