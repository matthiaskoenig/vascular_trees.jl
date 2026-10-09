module VascularTrees

    using Reexport # only for common modules, to avoid reexporting everything from the submodules

    include("common/paths.jl");        @reexport using .Paths
    include("common/definitions.jl");  @reexport using .Definitions
    include("common/options.jl");      @reexport using .Options
    include("common/table_utils.jl");   @reexport using .Options
    
    # part that processes initial trees (.vtk files) and saves them in .arrow format
    include("trees_processing/read_tree_files.jl")
    include("trees_processing/helpers.jl")
    include("trees_processing/process_individual_trees.jl")
    include("trees_processing/process_terminal_nodes.jl")
    
    # part that runs simulations
    include("models/pharmacokinetic_models.jl")

    include("simulation/get_ode_parameters.jl")
    include("simulation/benchmarking.jl")
    include("simulation/helpers.jl")

end