"""
`Altimetry`

Ocean surface altimetry information.
"""
struct Altimetry
	z::Function
	min::Float64
	max::Float64
end