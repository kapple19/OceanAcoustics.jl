"""
`Altimetry`

Ocean surface altimetry information.
"""
struct Altimetry
	z::Function
	min::Float64
	max::Float64

	function Altimetry(z::Function, min::Float64, max::Float64)
		zFcn(r::Real) = z(r)
		new(zFcn, min, max)
	end
end