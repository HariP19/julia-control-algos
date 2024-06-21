#=
The following code contains the code for finding 
    linear approximation of a function at a given
    point x0
=#
using ForwardDiff
using LinearAlgebra

function f(x)
    x1, x2, x3 = x 
    return [x1*x2*x3; log(2 + cos(x1)) + x2^x1; x1*x3/(1+ x2^2)]
end

# linear approximation at x0
# f(x) = f(x0) - J_x0 * (x - x0)
function f_lin_approx(x0, x)
    J_x0 = ForwardDiff.jacobian(f, x0)
    # display(J_x0)

    return f(x0) - J_x0 * (x - x0)
end

function main()
    x0 = [pi; 1.0; 2.0]
    perturb = 0.001
    x = x0 .+ perturb

    y = f(x)
    y_lin = f_lin_approx(x0, x)

    # Compute the relative error for each component
    rel_error = abs.(y - y_lin) ./ abs.(y)
    norm_error = norm(rel_error)

    println("Error between f and f_lin_approx at x0")
    display(norm_error)
end

# Call main
main()