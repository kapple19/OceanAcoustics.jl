using OceanAcoustics

scn = ExampleScenarios.convergence()

sols = OceanAcoustics.propagation(scn)

##
using Plots

p = plot()
plot!.(sols, vars = (1, 2))
display(p)
