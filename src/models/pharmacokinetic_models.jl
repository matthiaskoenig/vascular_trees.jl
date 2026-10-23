"""
Module with pharmacokinetic ODE based models for the drug(s): transport of the
    substance(s) along the liver's vessel tree (arterial/portal vascular trees ->
    liver sinusoids (terminal part) -> venous/biliary trees).

ODEs for inflow trees (arterial and portal) and for terminal part differ from ODEs for
    outflow trees (venous and biliary), so they are written in separate functions.
Also, variables (initial values and differentials) and parts of parameters structure 
    differ between ODE functions for vascular trees (vectors) and for terminal part 
    (arrays).
"""
module PharmacokineticModels

using ..Definitions: 
    TreeDefinitions, 
    TerminalParameters, 
    VascularTreeParameters

using InteractiveUtils, 
    SparseArrays

include("interventions.jl")
using .Interventions: f_intervention

include("equations/ode_trees.jl")
include("equations/ode_terminal.jl")
export jf_dxdt!

include("jacobians.jl")
export jac!

end
