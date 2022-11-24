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
	z_interp = linear_interpolation(r, z)
	z_fcn(r) = z_interp_fcn(r)
	Altimetry(z_fcn, minimum(z), maximum(z))
end