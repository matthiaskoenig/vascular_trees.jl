module ode_example
    export p, x0, f_dxdt

    # parameter vector
    p = [
        1.0, 
        5.0
    ]

    # initial conditions
    x0 = [
        10.0, 
        10.0
    ]

    #Define the problem
    function f_dxdt(dx, x, p, t)
        dx[1] = - p[1] * x[1]
        dx[2] = - p[2] * x[2]
    end

end