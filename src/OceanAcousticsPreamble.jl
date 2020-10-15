"An abstract supertype for the acoustic properties modelled in this module."
abstract type OceanAcoustic end

"Concise printing of the propreties of the Acoustic type and its subtypes"
Base.show(io::IO, oac::OceanAcoustic) = print_properties_types(io, oac)

"Enables broadcasting on the Acoustic type and its subtypes"
Base.broadcastable(oac::OceanAcoustic) = Ref(oac)
