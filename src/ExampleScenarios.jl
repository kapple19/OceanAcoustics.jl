export OAC_EXAMPLE_NAMES
export run_example

module ExampleScenarios
using OceanAcoustics
export flat

function flat()
	# Environment
	ocn = Medium(1500)
	bty = Boundary(5e2, 1600)
	env = Environment(10e2, ocn, bty)

	# Scenario
	fan = Fan(π/4 * LinRange(-1, 1, 51))
	src = Source(Position(0, 2e2), Signal(50), fan)
	return sno = Scenario(env, src)
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
