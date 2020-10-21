# Dependencies
using ForwardDiff: derivative
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
Rodas4
using DiffEqBase: AbstractODESolution

# OceanAcoustic Abstract Type
abstract type OceanAcoustic <: Any end

function Base.show(io::IO, oac::OceanAcoustic)
	println(io, string(typeof(oac)), "(")
	for p ∈ propertynames(oac)
		println(io, " ", p, "::", typeof(getproperty(oac, p)))
	end
	print(io, ")")
end

Base.broadcastable(oac::OceanAcoustic) = Ref(oac)