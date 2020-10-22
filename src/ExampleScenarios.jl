export OAC_EXAMPLE_NAMES
export trace_example

module ExampleScenarios
using OceanAcoustics
export flat
export slopes
export parabolic
export upward
export convergence
export wavy
export channel

function flat()
	# Environment
	ocn = Medium(1500)
	bty = Boundary(5e2, 1600)
	env = Environment(10e2, ocn, bty)

	# Scenario
	fan = Fan(π/4 * LinRange(-1, 1, 51))
	src = Source(Position(0, 2e2), Signal(50), fan)
	sno = Scenario(env, src, "Flat Environment")
end

function slopes()
	# Environment
	R = 10e3
	Z = 2e3
	c(r, z) = 1500 - 100r/R + 100z/Z
	zBty(r) = Z - 500r/R
	zAti(r) = 100r/R

	ocn = Medium(c)
	bty = Boundary(zBty, 1600)
	ati = Boundary(zAti, 343)
	env = Environment(R, ocn, bty, ati)

	# Scenario
	r₀ = 0
	z₀ = Z/4

	fan = Fan(acos(c(r₀, z₀)/c(r₀, Z)) * (-2:0.2:2))
	src = Source(Position(r₀, z₀), Signal(50), fan)
	sno = Scenario(env, src, "Sloped Environment")
end

function parabolic()
	# Environment
	c = 250
	b = 2.5e5
	zBty(r) = 2e-3b * √(1 + r/c)
	R = 20e3
	Z = 5e3

	ocn = Medium(c)
	bty = Boundary(zBty, 1600)
	env = Environment(R, ocn, bty)

	# Scenario
	fan = Fan(LinRange(atan(5e3/2e3), atan(5e3/20e3), 30))
	src = Source(Position(0, 0), Signal(50), fan)
	sno = Scenario(env, src, "Parabolic Bathymetry")
end

function upward()
	# Environment
	R = 1e5
	Z = 5e3
	c(r, z) = 1500 + 100z/Z

	ocn = Medium(c)
	bty = Boundary(Z, 1600)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 0

	fan = Fan(acos(c(r₀, z₀)/c(r₀, Z)) * LinRange(0.1, 1.0, 10))
	src = Source(Position(r₀, z₀), Signal(100), fan)
	sno = Scenario(env, src, "Upward-Refracting Rays")
end

function convergence()
	# Environment
	c = [1520, 1500, 1515, 1495, 1545.]
	z = [0., 300., 1200., 2e3, 5000.]
	Z = z[end]
	R = 250e3

	ocn = Medium((z, c))
	bty = Boundary(Z, 1600)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 20
	θ_crit = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, Z) |> acosd |> deg2rad

	fan = Fan(θ_crit * LinRange(0.5, 1, 11))
	src = Source(Position(r₀, z₀), Signal(200), fan)
	sno = Scenario(env, src, "Convergence-Zone Propagation")
end

# function upward_to_downward()
# 	# Environment
# 	R = 1e3
# 	Z = 5e3
# 	cr = [0, R]
# 	cz = [0, Z]
# 	c = []
# end

function wavy()
	# Altimetry
	zAtiMin = 0.
	zAtiMax = 50.
	zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2

	# Bathymetry
	rBtyPeak = 5e3
	zBtyMax = 1e3
	zBtyMin = 8e2
	Aᵣ = (2rBtyPeak/3)^2/log((9.0zBtyMax - 11.0zBtyMin)/(10.0(zBtyMax - zBtyMin)))
	zBty(r) = zBtyMax - (zBtyMax - zBtyMin)*exp(-(r - rBtyPeak)^2/4e5)

	# Ocean
	rOcnMax = 10e3
	cOcnMin = 1500.
	cOcnMax = 1600.
	cSolve(r) = [1.0 zAti(r) zAti(r)^2
		1.0 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2.0)^2
		1.0 zBty(r) zBty(r)^2]
	cSolved(r) = cSolve(r)\[cOcnMax, cOcnMin, cOcnMax]
	cCoeff₀(r) = cSolved(r)[1]
	cCoeff₁(r) = cSolved(r)[2]
	cCoeff₂(r) = cSolved(r)[3]
	cOcn(r, z) = cCoeff₀(r) + cCoeff₁(r)*z + cCoeff₂(r)*z^2

	# Environment
	ocn = Medium(cOcn)
	bty = Boundary(zBty, 1600)
	ati = Boundary(zAti, 343)
	env = Environment(rOcnMax, ocn, bty, ati)

	# Scenario
	r₀ = 0
	z₀ =(zBty(r₀) + zAti(r₀))/2

	fan = Fan(acos(cOcn(r₀, z₀)/cOcnMax).*(-1.5:0.125:1.5))
	src = Source(Position(r₀, z₀), Signal(200), fan)
	sno = Scenario(env, src, "Wavy Environment")
end

function channel()
	# Environment
	z = [0.0, 500/3, 500/2, 500, 1000, 1500, 4e3]
	c = [1480, 1500, 1485, 1475, 1480, 1485, 1525.]
	R = 250e3
	Z = z[end]

	ocn = Medium((z, c))
	bty = Boundary(Z, 1600)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 500
	θ_crit = ocn.SSP.c(0.0, z₀)/1500.0 |> acos

	fan = Fan(θ_crit * LinRange(-1, 1, 31))
	src = Source(Position(r₀, z₀), Signal(150), fan)
	sno = Scenario(env, src, "Deep Sound Channel")
end

end # ExampleScenarios

function example_names()
	example_names = names(OceanAcoustics.ExampleScenarios)
	filter!(name -> name ≠ :ExampleScenarios, example_names)
end

const OAC_EXAMPLE_NAMES = example_names()

function trace_example(name::Symbol)
	ex = getfield(ExampleScenarios, name)
	sno = ex()
	return trc = Trace(sno)
end
