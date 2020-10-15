"""
OceanAcoustics

A Julia package implementation of acoustic propagation in an ocean environment.
"""
module OceanAcoustics

"An abstract supertype for the acoustic properties modelled in this module."
abstract type OceanAcoustic end

"Concise printing of the propreties of the Acoustic type and its subtypes"
Base.show(io::IO, ac::OceanAcoustic) = print_properties_types(io, ac)

"Enables broadcasting on the Acoustic type and its subtypes"
Base.broadcastable(ac::OceanAcoustic) = Ref(ac)

include("AcousticPropagation.jl")

end