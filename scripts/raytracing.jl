## Preamble
using DrWatson
using OceanAcoustics
using Plots

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
		plot!(rays[nRay].Sol, vars = (1, 2))
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

## Flat Scenario
pars_rays_flat = Dict(
	:title => "Flat Environment",
	:z₀ => 200,
	:f => 200,
	:c => 1500,
	:R => 2000,
	:zBty => 500,
	:zAti => 0,
	:θ₀ => π/4*range(-1, 1, length = 51))

dict_rays_flat = Dict(
	:scen => "Flat",
	:pars => pars_rays_flat
)

## Smooth Scenario
# Altimetry
zAtiMin = -10
zAtiMax = 50
zAtiFcn(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2
# Bathymetry
rBtyPeak = 5e3
zBtyMax = 1e3
zBtyMin = 8e2
Aᵣ = (2rBtyPeak/3)^2/log((9zBtyMax - 11zBtyMin)/(10(zBtyMax - zBtyMin)))
zBtyFcn(r) = zBtyMax - (zBtyMax - zBtyMin)*exp(-(r - rBtyPeak)^2/4e5)
# Ocean
rOcnMax = 10e3
cOcnMin = 1500
cOcnMax = 1600
cSolve(r) = [1 zAtiFcn(r) zAtiFcn(r)^2
	1 (zAtiFcn(r) + zBtyFcn(r))/2 ((zAtiFcn(r) + zBtyFcn(r))/2)^2
	1 zBtyFcn(r) zBtyFcn(r)^2]
cSolved(r) = cSolve(r)\[cOcnMax, cOcnMin, cOcnMax]
cCoeff₀(r) = cSolved(r)[1]
cCoeff₁(r) = cSolved(r)[2]
cCoeff₂(r) = cSolved(r)[3]
cOcnFcn(r, z) = cCoeff₀(r) + cCoeff₁(r)*z + cCoeff₂(r)*z^2
# Source
r₀ = 0.
z₀ = (zBtyFcn(r₀) + zAtiFcn(r₀))/2
f = 250
# Rays
θ₀ = acos(cOcnFcn(r₀, z₀)/cOcnMax).*(-1.5:0.125:1.5)

pars_rays_smooth = Dict(
	:title => "Smooth Environment",
	:z₀ => z₀,
	:f => f,
	:c => cOcnFcn,
	:R => rOcnMax,
	:zBty => zBtyFcn,
	:zAti => zAtiFcn,
	:θ₀ => θ₀
)

dict_rays_smooth = Dict(
	:scen => "Smooth",
	:pars => pars_rays_smooth
)

## Parabolic Bathymetry
cOcnVal = 250
zBtyFcn(r) = 2e-3*2.5e5sqrt(1 + r/cOcnVal)
θ₀ = range(atan(5e3/2e3), atan(5e3/20e3), length = 30)

pars_rays_parabolic = Dict(
	:title => "Parabolic Bathymetry",
	:z₀ => 0,
	:f => 2e2,
	:c => cOcnVal,
	:R => 20e3,
	:zBty => zBtyFcn,
	:zAti => 0.0,
	:θ₀ => θ₀
)

dict_rays_parabolic = Dict(
	:scen => "Parabolic",
	:pars => pars_rays_parabolic
)

## Uniformly Increasing Celerity
cOcnFcn(r, z) = 1500 + 100z/5e3
θ_crit = acos(cOcnFcn(0, 0)/cOcnFcn(0, 5e3))

pars_rays_upward = Dict(
	:title => "Upward Refracting",
	:z₀ => 0,
	:f => 2e2,
	:c => (r, z) -> cOcnFcn,
	:R => 1e5,
	:zBty => 5e3,
	:zAti => 0.0,
	:θ₀ => θ_crit*range(0.1, 1, length = 10)
)

dict_rays_upward = Dict(
	:scen => "Upward",
	:pars => pars_rays_upward
)

## n²-Linear Profile
c₀ = 1550
z₀ = 1e3
cOcnFcn(r, z) = c₀/sqrt(1 + 2.4z/c₀)

pars_rays_n2linear = Dict(
	:title => "n²-Linear Profile",
	:z₀ => z₀,
	:f => 2e3,
	:c => cOcnFcn,
	:R => 3.5e3,
	:zBty => z₀,
	:zAti => 0.0,
	:θ₀ => acos(cOcnFcn(0, z₀)/cOcnFcn(0, 0))*range(-1.2, -0.4, length = 101)
)

dict_rays_n2linear = Dict(
	:scen => "n2Linear",
	:pars => pars_rays_n2linear
)

## Run Scenarios
# scenarios_rays = [
# 	dict_rays_flat,
# 	dict_rays_smooth,
# 	dict_rays_parabolic,
# 	dict_rays_upward,
# 	dict_rays_n2linear
# ]
scenarios_rays = [
	dict_rays_n2linear
]

run_sims(sims_rays, scenarios_rays)