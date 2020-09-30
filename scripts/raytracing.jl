## Preamble
using DrWatson
using OceanAcoustics
using Plots

## Load Scenarios
include("scenarios.jl")

## Scenario-Running Functions
function sim_rays(pars::Dict)
	@unpack z₀, f, c, R, zBty, zAti, θ₀, title = pars

	src = Source(Position(0, z₀), Signal(f))
	ocn = Medium(c, R)
	bty = Boundary(zBty)
	ati = Boundary(zAti)
	rays = Ray.(θ₀, src, ocn, bty, ati)

	pt = plot(
		xaxis = "Range (m)",
		yaxis = ("Depth (m)", :flip),
		title = "Ray Trace: " * title * " Scenario",
		legend = false
	)
	r = range(0, R, length = 1001)
	plot!(r, ati.z)
	plot!(r, bty.z)
	for nRay = 1:length(θ₀)
		plot!(rays[nRay].sol, vars = (1, 2))
	end
	return pt
end

function sims_rays(dict_rays::Dict)
	@unpack pars = dict_rays
	pt = sim_rays(pars)
	display(pt)
	return pt
end

function run_sims(sims_fcn::Function, scenarios)
	for (i, d) in enumerate(scenarios)
		p = sims_fcn(d)
		@show savename(d, "png")
		@show plotsdir("rays", savename(d, "png"))
		wsave(plotsdir("rays", savename(d, "png")), p)
	end
end

## Run Scenarios
scenarios_rays = [
	fcn_rays_flat(),
	fcn_rays_smooth(),
	fcn_rays_parabolic(),
	fcn_rays_upward(),
	fcn_rays_n2linear()
]

run_sims(sims_rays, scenarios_rays)