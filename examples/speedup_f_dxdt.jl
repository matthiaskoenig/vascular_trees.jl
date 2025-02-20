
# - use type annotations
# - use ! indicator on function
# - remove whitespace, comments and lines; 
# - vectorize calculation (matrix multiplications)

# tools:
# SnoopCompile: https://timholy.github.io/SnoopCompile.jl/stable/
# PreCompile: https://github.com/JuliaLang/PrecompileTools.jl
# PackageCompiler: https://github.com/JuliaLang/PackageCompiler.jl

# discussion: https://discourse.julialang.org/t/why-julia-compilation-times-are-so-long-in-case-of-functions-created-using-symbolics/112860
# discussion: Compile time and memory usage grows quickly with function size: https://github.com/JuliaLang/julia/issues/19158

# tips: Compiling functions with n = 40000 lines is may be a bad idea anyway from a performance perspective. At some point your function gets too big to fit in the instruction cache and then execution (not compilation) will slow down dramatically.
#       When I've seen others generate long function bodies and encounter this kind of problem, the solution was to automatically generate smaller functions with layers of calling functions over them.
# converting to c functions: https://symbolics.juliasymbolics.org/stable/tutorials/converting_to_C/; 
# see https://discourse.julialang.org/t/why-julia-compilation-times-are-so-long-in-case-of-functions-created-using-symbolics/112860/12

# I’ve seen people successfully abusing GCC to compile megabyte-sized math formulas in C. If it’s possible to translate your code to C, I’d try making a shared library with GCC and loading it in Julia. 
# https://symbolics.juliasymbolics.org/stable/manual/build_function/ : build_function with sharded_form
println("Starting include")
@time begin
    # @time include("Rectangle_trio_50_vectorized.jl")
    @time include("Rectangle_trio_50_symbolic.jl")
    @time @time_imports using .Transport_model

    # @time begin
    #     include("Rectangle_trio_1500_vectorized.jl")   # Load the odes module
    #     using .Transport_model         # Make the odes module available for use
    # end

    # println("First call")

    @time @eval Transport_model.f_dxdt!(
        zeros(size(Transport_model.x0)...),
        Transport_model.x0,
        Transport_model.p,
        0.0,
    )

    # println("Second call")
    @time @eval Transport_model.f_dxdt!(
        zeros(size(Transport_model.x0)...),
        Transport_model.x0,
        Transport_model.p,
        0.0,
    )

end;
