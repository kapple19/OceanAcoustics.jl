export ExampleScenarios
export OAC_EXAMPLE_NAMES
export example_scenario
export example_trace
export example_field
export example_grid

include("ExampleScenarios.jl")

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
	scn = example_scenario(name)
	trc = Trace(scn)
end

example_trace(name::String) = Symbol(name) |> example_trace

function example_field(name::Symbol)
	scn = example_scenario(name)
	fld = Field(scn)
end

example_field(name::String) = Symbol(name) |> example_field

function example_grid(name::Symbol)
	scn = example_scenario(name)
	grid = Grid(scn)
end

example_grid(name::String) = Symbol(name) |> example_grid
