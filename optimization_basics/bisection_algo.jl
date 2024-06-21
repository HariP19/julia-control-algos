function f(x)
    return (x-1)^2 - 4
end

# Bisection Params
global a = -1
global b = 5
global delta = 1e-9
global MAX_ITERS = 100
global N = 1

while N < MAX_ITERS
    global c = (a + b)/2; fc = f(c)
    if (abs(fc) < delta) || abs(c-a) < delta 
        print("Convergence reached at N = ", N, "\n")
        break
    end

    # Update for next iter 
    if sign(fc) == sign(f(a))
        global a = c 
    else
        global b = c
    end
    global N += 1
end

c = (a + b)/2
print("c: ", c, "\n")
print("f(c): ", f(c), "\n")