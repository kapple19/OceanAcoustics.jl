## Preamble
using OceanAcoustics
using DrWatson
using Plots

## Scenario-Running Functions
function sim_rays(scen::Function)
	θ₀, src, ocn, bty, ati, title = scen()

	rays = Ray.(θ₀, src, ocn, bty, ati)
	r = range(0, ocn.R, length = 1001)

	pt = plot(
		xaxis = "Range (m)",
		yaxis = ("Depth (m)", :flip),
		title = "Ray Trace: " * title,
		legend = false
	)
	plot!(r, ati.z)
	plot!(r, bty.z)
	for nRay = 1:length(θ₀)
		plot!(rays[nRay].sol, vars = (1, 2))
	end
	return pt
end

# function sims_rays(dict_rays::Vector)
# 	@unpack pars = dict_rays
# 	pt = sim_rays(pars)
# 	display(pt)
# 	return pt
# end

function run_sims(sim_fcn::Function, scenarios::AbstractVector)
	for scen ∈ scenarios
		p = sim_fcn(scen)
		filename = String(Symbol(scen))
		wsave(plotsdir("rays") * "/" * filename * ".png", p)
	end
end

## Run Scenarios
include("scenarios.jl")

# scenarios_rays = [
# 	flat,
# 	smooth,
#	convergence,
#	upward
# ]

scenarios_rays = [upward]

run_sims(sim_rays, scenarios_rays)