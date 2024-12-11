using OrdinaryDiffEq, Plots

include("odes.jl")   # Load the odes module
using .ode_example          # Make the odes module available for use

# Ode integration
tspan = (0.0, 10.0)
prob = ODEProblem(ode_example.f_dxdt, ode_example.x0, tspan, ode_example.p)
sol = solve(prob, Tsit5())

print(sol)
