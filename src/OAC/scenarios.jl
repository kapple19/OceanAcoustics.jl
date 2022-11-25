"""
`Depth`

Storing univariate real-valued function (F64 -> F64) and meta-information.

`Depth(z::Function, min::Real, max::Real)`
`Depth(r::Vector{<:Real}, z::Vector{<:Real})` creates an interpolator
`Depth(z::Real)` creates a function

Used for:
	* ocean surface altimetry
	* ocean bottom bathymetry
"""
mutable struct Depth
	fcn::Function
	min::Float64
	max::Float64

	Depth(z::Function) = new(z)
end

export Depth

function Depth(fcn::Function, min::Float64, max::Float64)
	d = Depth(fcn)
	d.min = min
	d.max = max
	return d
end

function Depth(fcn::Function, min::Real, max::Real)
	zFcn(r::Real) = fcn(r)
	Depth(zFcn, Float64(min), Float64(max))
end

function Depth(r::Vector{<:Real}, z::Vector{<:Real})
	# Checks
	!issorted(r) && throw(NotSorted(r))
	!allunique(r) && throw(NotAllUnique(r))
	length(r) ≠ length(z) && throw(DimensionMismatch())

	Depth(linear_interp_fcn(r, z), minimum(z), maximum(z))
end

function Depth(z::Real)
	zF64 = Float64(z)
	Depth(r -> zF64, zF64, zF64)
end

(z::Depth)(r) = z.fcn(r)