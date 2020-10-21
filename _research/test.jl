## Simple Example
using OceanAcoustics
using Plots

# Environment
ocn = Medium(1500)
bty = Boundary(1e3, 1600)
env = Environment(5e3, ocn, bty)

# Scenario
spk = Spark(π/4)
fan = Fan(spk)
src = Source(Position(0, 2e2), Signal(50), fan)
sno = Scenario(env, src)

# Trace
trc = Trace(sno)

# Plot
pt = plot(legend = false)
plot!.([trc.rays[nRay].sol for nRay = eachindex(trc.rays)], vars = (1, 2))
display(pt)

## Running built-in examples
using OceanAcoustics
using Plots

trcs = run_example.(OAC_EXAMPLE_NAMES)

for nTrc ∈ eachindex(trcs)
	pt = plot(legend = false)
	plot!.([trcs[nTrc].rays[nRay].sol for nRay = eachindex(trcs[nTrc].rays)], vars = (1, 2), yaxis = :flip)
	display(pt)
end
