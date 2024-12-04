using OrdinaryDiffEq
using Plots
using DataFrames
using CSV
using BenchmarkTools
using Distributed
using Sundials

include("Rectangle_quad_10.jl")   # Load the odes module
using .Rectangle_quad_10          # Make the odes module available for use

println("Number of threads: ", Threads.nthreads())
println("Number of processes: ", nprocs())

absolute_tolerance = 1e-6
relative_tolerance = 1e-6

# x0[26] ged_A_marginal
x0 = Rectangle_quad_10.x0
x0[26] = 1.0  # set marginal concentration

# Ode integration
tspan = (0.0, 10.0/60)
tpoints = range(tspan[1], stop=tspan[2], length=1001)
prob = ODEProblem(
    Rectangle_quad_10.f_dxdt, 
    x0,
    tspan, 
    Rectangle_quad_10.p)

@time begin
sol = solve(
    prob, 
    CVODE_BDF(), # Rosenbrock23(), # Tsit5(), # CVODE_BDF
    saveat=tpoints,
    dense=false,
    reltol=relative_tolerance, 
    abstol=absolute_tolerance)
end

# print(sol)

# Step 3: Convert solution to DataFrame
df = DataFrame(time=sol.t, value=sol.u)

# Step 4: Write DataFrame to CSV
header = vcat(["time"], Rectangle_quad_10.xids)
CSV.write("Rectangle_quad_10.csv", df, header=header)
