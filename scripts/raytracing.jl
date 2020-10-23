## Preamble
using OceanAcoustics
using DrWatson
using GRUtils

## Simulate Rays
function sim_rays(scen::Function)
	θ₀, src, ocn, bty, ati, title = scen()
	println("Simulating Rays for ", title)

	rays = Ray.(θ₀, src, ocn, bty, ati)
	
	f = acoustic_plot(ati)
	acoustic_plot!(bty)
	acoustic_plot!.(rays)
	acoustic_plot!("Ray Trace: " * title)
	return f
end

## Run Scenarios
include("scenarios.jl")
include("simulations.jl")

# scenarios_rays = [
# 	flat,
# 	smooth,
# 	convergence,
# 	upward,
# 	parabolic,
# 	channel,
# 	seamount,
# 	simple,
# 	n2linear,
# 	slopes
# ]

scenarios_rays = [
	smooth,
	convergence
]


run_sims(sim_rays, scenarios_rays)
