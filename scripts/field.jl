## Preamble
using OceanAcoustics
using DrWatson
using Plots

## Simulate Field
function sim_field(scen::Function)
	θ₀, src, ocn, bty, ati, title = scen()

	fld = Field(θ₀, src, ocn, bty, ati)

	rng = range(0, ocn.R, length = 51)
	dpt = range(0, ocn.Z, length = 31)

	p = heatmap(rng, dpt, fld.TL,
		xaxis = "Range (m)",
		yaxis = ("Depth (m)", :flip))
	return p
end

## Run Scenarios
include("scenarios.jl")
include("simulations.jl")

# scenarios_field = [
# 	convergence
# ]

scenarios_field = [
	n2linear
]

run_sims(sim_field, scenarios_field)