# Testing the import of SBML models for simulations

using SBMLToolkit, OrdinaryDiffEq
using Sundials

# odesys = readSBML("src/sbml/dapagliflozin_body_flat.xml", ODESystemImporter())
odesys = readSBML("src/sbml/dapagliflozin_body_flat.xml", ReactionSystemImporter())
# odesys = readSBML("src/sbml/dapagliflozin_liver.xml", ODESystemImporter())
# odesys = readSBML("src/sbml/dapagliflozin_liver.xml", ReactionSystemImporter())


tspan = (0.0, 24 * 60.0)
tpoints = range(tspan[1], stop=tspan[2], length=1001)

prob = ODEProblem(odesys, [], tspan, [])
@time begin
sol = solve(
    prob, 
    # Tsit5(), 
    CVODE_BDF(), # Rosenbrock23(), # Tsit5(), # CVODE_BDF
    dense=false, saveat=tpoints)
end
