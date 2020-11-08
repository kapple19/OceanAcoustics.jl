using ForwardDiff:
derivative,
gradient
using Interpolations:
LinearInterpolation,
Flat
using IntervalArithmetic:
Interval,
(..)
using OrdinaryDiffEq:
ODEProblem,
solve,
ContinuousCallback,
CallbackSet,
terminate!,
AutoVern7,
Rodas4,
ODECompositeSolution
using Roots: find_zeros
using ProgressMeter: @showprogress

"""
`OceanAcoustic <: Any`

An abstract supertype for all `struct`s defined for the ocean acoustics modelling implementation.
"""
abstract type OceanAcoustic <: Any end

function Base.show(io::IO, oac::OceanAcoustic)
	println(io, string(typeof(oac)), "(")
	for p ∈ propertynames(oac)
		println(io, " ", p, "::", typeof(getproperty(oac, p)))
	end
	print(io, ")")
end

Base.broadcastable(oac::OceanAcoustic) = Ref(oac)
