

"""
`Depth`

Storing univariate real-valued function (F64 -> F64).

`Depth(z::Function)`
`Depth(z::Function, min::Real, max::Real)`
`Depth(r::Vector{<:Real}, z::Vector{<:Real})` creates an interpolator
`Depth(z::Real)` creates a function

Used for:
	* ocean surface altimetry
	* ocean bottom bathymetry
"""
struct Depth <: Oac
	fcn::Function

	function Depth(z::Function)
		fcn(r::Real) = Float64(z(r))
		new(fcn)
	end
end

(z::Depth)(r) = z.fcn(r)

function Depth(r::Vector{<:Real}, z::Vector{<:Real})
	# Checks
	!issorted(r) && throw(NotSorted(r))
	!allunique(r) && throw(NotAllUnique(r))
	length(r) ≠ length(z) && throw(DimensionMismatch())

	z_interp = linear_interpolation(r, z, extrapolation_bc = Line())
	z_fcn(r::Real) = z_interp(r)

	Depth(z_fcn)
end

function Depth(z::Real)
	fcn(r::Real) = Float64(z)
	Depth(fcn)
end

Depth(args) = Depth(args...)

# function effective_acoustic_admittance(ρ::Real, c::Number, θ)
# 	A = sin(θ) / ρ / c
# end

# function surface_reflection_coefficient(
# 	ρ₁::Number,
# 	ρ₂::Number,
# 	c₁::Number,
# 	c₂::Number,
# 	θ₁::Number,
# 	θ₂::Number
# )
# 	A₁ = effective_acoustic_admittance(ρ₁, c₁, θ₁)
# 	A₂ = effective_acoustic_admittance(ρ₂, c₂, θ₂)
# 	R = (A₁ - A₂) / (A₁ + A₂)
# end

# function bottom_reflection_coefficient(
# 	ρ₁::Number,
# 	ρ₂::Number,
# 	c₁::Number,
# 	cₚ::Number,
# 	cₛ::Number,
# 	θ₁::Number
# )

# 	θₚ = acos(cₚ / c₁ * cos(θ₁))
# 	θₛ = acos(cₛ / c₁ * cos(θ₁))
	
# 	Z₁ = 1 / effective_acoustic_admittance(ρ₁, c₁, θ₁)
# 	Zₚ = 1 / effective_acoustic_admittance(ρₚ, cₚ, θₚ)
# 	Zₛ = 1 / effective_acoustic_admittance(ρₛ, cₛ, θₛ)
# 	Zₜ = Zₚ * cos(2θₛ)^2 + Zₛ * sin(2θₛ)^2
# 	R = (Zₜ - Z₁) / (Zₜ + Z₁)
# end

# """
# `ReflectionCoefficient`
# """
# struct ReflectionCoefficient <: Oac
# 	fcn::Function
# end

# (R::ReflectionCoefficient)(θ) = R.fcn(θ)

# function ReflectionCoefficient(R::Number)
# 	fcn(θ::Number) = R
# end

# """
# Surface `ReflectionCoefficient`
# """
# function ReflectionCoefficient(
# 	ρ₁::Number,
# 	ρ₂::Number,
# 	c₁::Number,
# 	cₚ::Number,
# 	cₛ::Number,
# 	θ₁::Number
# )
# 	R = surface_reflection_coefficient(ρ₁, ρ₂, c₁, cₚ, cₛ, θ₁)
# 	ReflectionCoefficient(R)
# end

# """
# Bottom `ReflectionCoefficient`
# """
# function ReflectionCoefficient(
# 	ρ₁::Number,
# 	ρ₂::Number,
# 	c₁::Number,
# 	cₚ::Number,
# 	cₛ::Number,
# 	θ₁::Number
# )
# 	R = bottom_reflection_coefficient(ρ₁, ρ₂, c₁, cₚ, cₛ, θ₁)
# 	ReflectionCoefficient(R)
# end

# ReflectionCoefficient(args) = ReflectionCoefficient(args...)

# struct Celerity <: Oac
# 	fcn::Function
# end

# (c::Celerity)(r, z) = c.fcn(r, z)

