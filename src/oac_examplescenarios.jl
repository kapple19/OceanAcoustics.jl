export ExampleScenarios
export OAC_EXAMPLE_NAMES
export example_scenario
export example_trace
export example_field
export example_grid

module ExampleScenarios
using OceanAcoustics

function convergence()
	c = [1520, 1500, 1515, 1495, 1545.]
	z = [0., 300., 1200., 2e3, 5000.]
	Z = z[end]
	R = 250e3

	ocn = Medium(z, c)
	bty = Boundary(5e3)
	ati = Boundary(0.)
	env = Environment(R, ocn, bty, ati)

	θ_crit = acos(ocn.SSP.c(0.0, 0.0)/ocn.SSP.c(0.0, 5e3))
	θ₀s = θ_crit*LinRange(0.5, 1.0, 10)

	fan = Fan(θ₀s)
	src = Source(Position(0., 0.), Signal(200.), fan)
	scn = Scenario(env, src, "Convergence Zone Propagation")
end

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
	cSolve(r) = [
		1.0 zAti(r) zAti(r)^2
		1.0 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2)^2
		1.0 zBty(r) zBty(r)^2
	]
	cCoeff(r) = cSolve(r)\[cOcnMax, cOcnMin, cOcnMax]
	cOcn(r, z) = cCoeff(r)' * [1, z, z^2]

	# Environment
	ocn = Medium(cOcn)
	bty = Boundary(zBty)
	ati = Boundary(zAti)
	env = Environment(rOcnMax, ocn, bty, ati)

	# Scenario
	r₀ = 0.0
	z₀ = (zBty(r₀) + zAti(r₀))/2
	θ_crit = acos(cOcn(r₀, z₀)/cOcnMax)

	fan = Fan(θ_crit*(-1.5:0.125:1.5))
	src = Source(Position(r₀, z₀), Signal(250.), fan)
	scn = Scenario(env, src, "Wavy Environment")
end

function flat()
	# Environment
	ocn = Medium(1500)
	bty = Boundary(5e2)
	env = Environment(10e2, ocn, bty)

	# Scenario
	fan = Fan(π/4 * LinRange(-1, 1, 51))
	src = Source(Position(0, 2e2), Signal(50), fan)
	scn = Scenario(env, src, "Flat Environment")
end

function slopes()
	# Environment
	R = 10e3
	Z = 2e3
	c(r, z) = 1500 - 100r/R + 100z/Z
	zBty(r) = Z - 500r/R
	zAti(r) = 100r/R

	ocn = Medium(c)
	bty = Boundary(zBty)
	ati = Boundary(zAti)
	env = Environment(R, ocn, bty, ati)

	# Scenario
	r₀ = 0
	z₀ = Z/4

	fan = Fan(acos(c(r₀, z₀)/c(r₀, Z)) * (-2:0.2:2))
	src = Source(Position(r₀, z₀), Signal(50), fan)
	scn = Scenario(env, src, "Sloped Environment")
end

function parabolic()
	# Environment
	c = 250
	b = 2.5e5
	zBty(r) = 2e-3b * √(1 + r/c)
	R = 20e3
	Z = 5e3

	ocn = Medium(c)
	bty = Boundary(zBty)
	env = Environment(R, ocn, bty)

	# Scenario
	fan = Fan(LinRange(atan(5e3/2e3), atan(5e3/20e3), 30))
	src = Source(Position(0, 0), Signal(50), fan)
	scn = Scenario(env, src, "Parabolic Bathymetry")
end

function upward()
	# Environment
	R = 1e5
	Z = 5e3
	c(r, z) = 1500 + 100z/Z

	ocn = Medium(c)
	bty = Boundary(Z)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 0

	fan = Fan(acos(c(r₀, z₀)/c(r₀, Z)) * LinRange(0.1, 1.0, 10))
	src = Source(Position(r₀, z₀), Signal(100), fan)
	scn = Scenario(env, src, "Upward-Refracting Rays")
end

function seamount()
	# Environment
	zc = [0, 100, 200, 350, 500, 1500, 3100.]
	c = [1480, 1470, 1475, 1473, 1475, 1488, 1505.]

	Z = zc[end]

	rBty = 1e3*[0, 40, 45, 50, 55, 60, 70, 140]
	zBty = [Z, Z, 2900, 2850, 2000, 500, Z, Z]

	R = rBty[end]

	ocn = Medium(zc, c)
	bty = Boundary(rBty, zBty)
	env = Environment(R, ocn, bty)

	# Scenario
	fan = Fan(atan(363/2e3) * LinRange(-1, 1, 31))
	src = Source(Position(0, 363), Signal(200), fan)
	scn = Scenario(env, src, "Seamount")
end

function channel()
	# Environment
	z = [0.0, 500/3, 500/2, 500, 1000, 1500, 4e3]
	c = [1480, 1500, 1485, 1475, 1480, 1485, 1525.]
	R = 250e3
	Z = z[end]

	ocn = Medium(z, c)
	bty = Boundary(Z)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 500
	θ_crit = ocn.SSP.c(r₀, z₀)/1500.0 |> acos

	fan = Fan(θ_crit * LinRange(-1, 1, 31))
	src = Source(Position(r₀, z₀), Signal(150), fan)
	scn = Scenario(env, src, "Deep Sound Channel")
end

function munk()
	# Environment
	z_(z) = 2(z - 1300)/1300
	ϵ = 7.37e-3
	c(r, z) = 1500*(1 + ϵ*(z_(z) - 1 + exp(-z_(z))))

	ocn = Medium(c)
	bty = Boundary(5e3)
	env = Environment(100e3, ocn, bty)

	# Scenario
	r₀ = 0.0
	z₀ = 1e3
	θ_crit = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, 0) |> acos

	fan = Fan(θ_crit * LinRange(-1, 1, 31))
	src = Source(Position(r₀, z₀), Signal(150), fan)
	scn = Scenario(env, src, "Munk Profile")
end

function n2linear()
	# Environment
	c₀ = 1550
	c(r, z) = c₀/√(1 + 2.4z/c₀)
	R = 4e3
	Z = 1e3

	ocn = Medium(c)
	bty = Boundary(Z)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0.0
	z₀ = Z
	θ_crit = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, 0) |> acos

	fan = Fan(θ_crit * LinRange(-1.2, -0.8, 21))
	src = Source(Position(r₀, z₀), Signal(2e3), fan)
	scn = Scenario(env, src, "n²-Linear Profile")
end

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end # Module: ExampleScenarios

function example_names()
	example_names = names(OceanAcoustics.ExampleScenarios)
	filter!(name -> name ≠ :ExampleScenarios, example_names)
end

const OAC_EXAMPLE_NAMES = example_names()

function example_scenario(name::Symbol)
	xmp = getfield(ExampleScenarios, name)
	scn = xmp()
end

example_scenario(name::String) = Symbol(name) |> example_scenario

function example_trace(name::Symbol)
	println("Ray Trace: " * String(name))
	scn = example_scenario(name)
	trc = Trace(scn)
end

example_trace(name::String) = Symbol(name) |> example_trace

function example_field(name::Symbol)
	println("Field: " * String(name))
	scn = example_scenario(name)
	fld = Field(scn)
end

example_field(name::String) = Symbol(name) |> example_field

function example_grid(name::Symbol)
	println("Grid: " * String(name))
	scn = example_scenario(name)
	grid = Grid(scn)
end

example_grid(name::String) = Symbol(name) |> example_grid
