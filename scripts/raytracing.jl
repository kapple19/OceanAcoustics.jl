## Preamble
using OceanAcoustics
using DrWatson
using GRUtils

## Simulate Rays
function sim_rays(scen::Function)
	θ₀, src, ocn, bty, ati, title = scen()

	rays = Ray.(θ₀, src, ocn, bty, ati)
	rng = range(0, ocn.R, length = 1001)

	f = acoustic_plot(ati)
	acoustic_plot!(f, bty)
	acoustic_plot!.(f, rays)
	acoustic_plot!(f, title)
	display(f)
	return f
end

## Run Scenarios
include("scenarios.jl")
include("simulations.jl")

scenarios_rays = [
	flat,
	smooth,
	convergence,
	upward,
	parabolic,
	channel,
	seamount,
	simple,
	n2linear
]

# scenarios_rays = [
# 	parabolic
# ]

run_sims(sim_rays, scenarios_rays)
