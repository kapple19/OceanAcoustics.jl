"""
OceanAcoustics

A Julia package implementation of acoustic propagation in an ocean environment.
"""
module OceanAcoustics

include("Preamble.jl")
include("MathSupport.jl")
include("Propagation.jl")
include("OACPlots.jl")
include("ExampleScenarios.jl")

end
