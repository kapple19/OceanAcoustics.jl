"""
OceanAcoustics

A Julia package implementation of acoustic propagation in an ocean environment.
"""
module OceanAcoustics

# export ExampleScenarios

include("OceanAcousticsPreamble.jl")
include("AugmentingAdmin.jl")
include("AugmentingMaths.jl")
include("AcousticPropagation.jl")
# include("ExampleScenarios.jl")

end