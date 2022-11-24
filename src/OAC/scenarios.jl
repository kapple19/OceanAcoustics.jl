"""
`Altimetry`

Ocean surface altimetry information.
"""
struct Altimetry
	z::Function
	min::Float64
	max::Float64

	function Altimetry(z::Function, min::Real, max::Real)
		zFcn(r::Real) = z(r)
		new(zFcn, Float64(min), Float64(max))
	end
end

export Altimetry

function Altimetry(r::Vector{<:Real}, z::Vector{<:Real})
	!issorted(r) && throw(NotSorted(r))
	!allunique(r) && throw(NotAllUnique(r))
	length(r) ≠ length(z) && throw(DimensionMismatch())
	Altimetry(linear_interp_fcn(r, z), minimum(z), maximum(z))
end

function Altimetry(z::Real)
	zF64 = Float64(z)
	Altimetry(r -> zF64, zF64, zF64)
end