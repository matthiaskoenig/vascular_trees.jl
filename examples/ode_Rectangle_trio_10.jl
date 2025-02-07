using OrdinaryDiffEq
using DataFrames
using CSV
using BenchmarkTools
using Distributed
using Sundials
using ParameterizedFunctions

@time begin
    include("Rectangle_trio_10.jl")   # Load the odes module
    using .Rectangle_trio_10          # Make the odes module available for use
end

@time begin
    include("Rectangle_trio_10_symbolic.jl")   # Load the odes module
    using .Transport_model         # Make the odes module available for use
end

@time begin
    include("Rectangle_trio_10_vectorized.jl")   # Load the odes module
    using .Transport_model         # Make the odes module available for use
end

@time begin
    include("Rectangle_trio_50_vectorized.jl")   # Load the odes module
    using .Transport_model         # Make the odes module available for use
end

@time begin
    include("Rectangle_trio_1500_vectorized.jl")   # Load the odes module
    using .Transport_model         # Make the odes module available for use
end


println("Number of threads: ", Threads.nthreads())
println("Number of processes: ", nprocs())

absolute_tolerance = 1e-6
relative_tolerance = 1e-6

# x0[26] ged_A_marginal
x0 = Transport_model.x0
x0[26] = 1.0  # set marginal concentration

# Ode integration
tspan = (0.0, 10.0 / 60)
tpoints = range(tspan[1], stop = tspan[2], length = 1001)

dx = zeros(size(x0)...)
@time begin
    Transport_model.f_dxdt(
        zeros(size(Transport_model.x0)...),
        Transport_model.x0,
        Transport_model.p,
        0.0,
    )
end

@time begin
    prob = ODEProblem(Transport_model.f_dxdt, Transport_model.x0, tspan, Transport_model.p)
end

@time begin
    sol = solve(
        prob,
        # CVODE_BDF(), # Rosenbrock23(), # Tsit5(), # CVODE_BDF
        Tsit5(),
        saveat = tpoints,
        dense = false,
        reltol = relative_tolerance,
        abstol = absolute_tolerance,
    )
end

x0[26] = 10.0  # set marginal concentration
@time begin
    sol = solve(
        prob,
        CVODE_BDF(), # Rosenbrock23(), # Tsit5(), # CVODE_BDF
        # Tsit5(),
        saveat = tpoints,
        dense = false,
        reltol = relative_tolerance,
        abstol = absolute_tolerance,
    )
end

# print(sol)

# Step 3: Convert solution to DataFrame
df = DataFrame(time = sol.t, value = sol.u)

# Step 4: Write DataFrame to CSV
header = vcat(["time"], Rectangle_trio_10.xids)
CSV.write("Rectangle_quad_10.csv", df, header = header)
