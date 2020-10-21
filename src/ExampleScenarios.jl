export OAC_EXAMPLE_NAMES
export run_example

module ExampleScenarios
using OceanAcoustics
export flat
export slopes

function flat()
	# Environment
	ocn = Medium(1500)
	bty = Boundary(5e2, 1600)
	env = Environment(10e2, ocn, bty)

	# Scenario
	fan = Fan(π/4 * LinRange(-1, 1, 51))
	src = Source(Position(0, 2e2), Signal(50), fan)
	sno = Scenario(env, src)
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
	sno = Scenario(env, src)
end

end # ExampleScenarios

function example_names()
	example_names = names(OceanAcoustics.ExampleScenarios)
	filter!(name -> name ≠ :ExampleScenarios, example_names)
end

const OAC_EXAMPLE_NAMES = example_names()

function run_example(name::Symbol)
	ex = getfield(ExampleScenarios, name)
	sno = ex()
	return trc = Trace(sno)
end
