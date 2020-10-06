## Preamble
using OceanAcoustics
using DrWatson
using Plots

## Simulate Field
function sim_field(scen::Function)
	θ₀, src, ocn, bty, ati, title = scen()
	rng = range(0, ocn.R, length = 51)
	dpt = range(0, ocn.Z, length = 31)

	fld = @time Field(θ₀, rng, dpt, src, ocn, bty, ati)

	TL = min.(100, -20log10.(abs.(fld.p)))'
	p = heatmap(rng, dpt, TL,
		xaxis = "Range (m)",
		yaxis = ("Depth (m)", :flip))
	return p
end

## Run Scenarios
include("scenarios.jl")
include("simulations.jl")

scenarios_field = [
	convergence
]

run_sims(sim_field, scenarios_field)