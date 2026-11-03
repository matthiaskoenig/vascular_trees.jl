using OrdinaryDiffEq, Plots
gr()

# parameter vector
p = [1.0, 5.0]

# initial conditions
x0 = [10.0, 10.0]

#Define the problem
function odes(dx, x, p, t)
    dx[1] = -p[1] * x[1]
    dx[2] = -p[2] * x[2]
end

# Ode integration
tspan = (0.0, 10.0)
prob = ODEProblem(odes, x0, tspan, p)
sol = solve(prob, Tsit5())

print(sol)

# # Plot
# plot(sol, linewidth = 2, title = "ODE Test",
#     xaxis = "Time", yaxis = "x",
#     label = "Numerical Solution"
# )
