## Preamble
using OceanAcoustics
using DrWatson
using Plots

## Simulate Field
function sim_field(scen::Function)
	θ₀, src, ocn, bty, ati, title = scen()

	fld = Field(θ₀, src, ocn, bty, ati)

	rng = range(0, ocn.R, length = 5)
	dpt = range(0, ocn.Z, length = 3)

	p = acoustic_plot(rng, dpt, fld)
	acoustic_plot!(rng, ati)
	acoustic_plot!(rng, bty)
	acoustic_plot!(extrema(rng), (0., ocn.Z))
	title!(title)
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