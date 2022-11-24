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

function Altimetry(r::Vector{<:Real}, z::Vector{<:Real})
	# check dimensions
	# check r is sorted
	Altimetry(linear_interp_fcn(r, z), minimum(z), maximum(z))
end