## Preamble
using OceanAcoustics
using DrWatson
using Plots

## Simulate Rays
function sim_rays(scen::Function)
	θ₀, src, ocn, bty, ati, title = scen()

	rays = Ray.(θ₀, src, ocn, bty, ati)
	rng = range(0, ocn.R, length = 1001)

	p = acoustic_plot(rays)
	acoustic_plot!(rng, ati)
	acoustic_plot!(rng, bty)
	acoustic_plot!(extrema(rng), (0., ocn.Z))
	title!(title)
	return p
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

run_sims(sim_rays, scenarios_rays)