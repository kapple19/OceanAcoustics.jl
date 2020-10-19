module OceanAcoustics
## Preamble
abstract type OceanAcoustic end

function Base.show(io::IO, oac::OceanAcoustic)
	println(io, string(typeof(oac)), "(")
	for p ∈ propertynames(oac)
		println(io, " ", p, "::", typeof(getproperty(oac, p)))
	end
	print(io, ")")
end

Base.broadcastable(oac::OceanAcoustic) = Ref(oac)

## Augmenting Mathematics
using ForwardDiff: derivative
using Interpolations:
LinearInterpolation,
Flat

function interpolated_function(x, y)
	itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
	return itp_fcn(x::Real) = itp(x)
end
function interpolated_function(x, y, z)
	itp = LinearInterpolation((x, y), z, extrapolation_bc = Flat())
	return itp_fcn(x::Real, y::Real) = itp(x, y)
end

function bivariate_derivatives(f::Function)
	∂f_∂x(x, y) = derivative(x -> f(x, y), x)
	∂f_∂y(x, y) = derivative(y -> f(x, y), y)

	∂²f_∂x²(x, y) = derivative(x -> ∂f_∂x(x, y), x)
	∂²f_∂y∂x(x, y) = derivative(y -> ∂f_∂x(x, y), y)
	∂²f_∂x∂y(x, y) = derivative(x -> ∂f_∂y(x, y), x)
	∂²f_∂y²(x, y) = derivative(y -> ∂f_∂y(x, y), y)

	return ∂f_∂x, ∂f_∂y, ∂²f_∂x², ∂²f_∂x∂y, ∂²f_∂y∂x, ∂²f_∂y²
end

## Acoustic Propagation
using ForwardDiff: derivative
using IntervalArithmetic: Interval, (..)

export (..)
export Boundary
export Celerity
export Medium
export Environment

function boundary_reflection(t_inc::Vector, t_bnd::Vector)
	MyAngle(tng) = atan(tng[2]/tng[1])
	θ_inc = MyAngle(t_inc)
	θ_bnd = MyAngle(t_bnd)

	c = cos(θ_inc)/t_inc[1]

	θ_inc_flat = θ_inc - θ_bnd
	θ_rfl_flat = -θ_inc_flat
	θ_rfl = θ_rfl_flat + θ_bnd

	return [cos(θ_rfl), sin(θ_rfl)]/c
end

struct Boundary <: OceanAcoustic
	depth::Depth
	c::Function
	callback::ContinuousCallback

	function Boundary(z::Function, c::Function)
		dz_dr(r) = derivative(z, r)

		condition(u, t, ray) = z(u[1]) - u[2]
		function affect!(ray)
			rᵢ = ray.u[1]
			ξᵢ = ray.u[3]
			ζᵢ = ray.u[4]

			ξₒ, ζₒ = boundary_reflection(
				[ξᵢ, ζᵢ],
				[1, dz_dr(rᵢ)]
			)

			if ξₒ < 0
				return terminate!(ray)
			else
				function reflect!(ray)
					ray.u[3] = ξₒ
					ray.u[4] = ζₒ
				end
				return reflect!(ray)
			end
		end
		callback = ContinuousCallback(condition, affect!)

		return new(z, c, callback)
	end
end

function Boundary(rz::AbstractArray{T}) where T <: Real
	r = [rng for rng ∈ rz[:, 1]]
	z = [dpt for dpt ∈ rz[:, 2]]
	return Boundary(r, z)
end

function Boundary(z::Real)
	zFcn(r) = z
	return Boundary(zFcn)
end

struct Celerity <: OceanAcoustic
	c::Function
	∂c_∂r::Function
	∂c_∂z::Function
	∂²c_∂r²::Function
	∂²c_∂r∂z::Function
	∂²c_∂z∂r::Function
	∂²c_∂z²::Function

	function Celerity(c::Function)
		∂c_∂r, ∂c_∂z, ∂²c_∂r², ∂²c_∂r∂z, ∂²c_∂z∂r, ∂²c_∂z² = bivariate_partial_derivatives(c)
		return new(
			c, ∂c_∂r, ∂c_∂z,
			∂²c_∂r²,
			∂²c_∂r∂z, ∂²c_∂z∂r,
			∂²c_∂z²)
	end
end

function Celerity(z::AbstractVector{T}, c::AbstractVector{T}) where T <: Real
	cFcn_z = interpolated_function(z, c)
	cFcn(r, z) = cFcn_z(z)
	return Celerity(cFcn)
end

function Celerity(c::Real)
	cFcn(r, z) = c
	return Celerity(cFcn)
end

struct Medium <: OceanAcoustic
	SSPₚ::Celerity
	SSPₛ::Celerity
end

function Medium(SSPₚ::Celerity)
	SSPₛ = Celerity(0)
	return Medium(SSPₚ, SSPₛ)
end

function Medium(c::Union{Function, Real})
	SSPₚ = Celerity(c)
	return Medium(SSPₚ)
end

struct Environment
	Ωr::Interval
	Ωz::Interval
	ocn::Medium
	bty::Boundary
	ati::Boundary

	function Environment(
		Ωr::Interval,
		ocn::Medium,
		bty::Boundary,
		ati::Boundary = Boundary(0))

		Ωz_ati = ati.z(Ωr)
		Ωz_bty = bty.z(Ωr)
		Ωz = Ωz_ati ∪ Ωz_bty

		return new(Ωr, Ωz, ocn, bty, ati)
	end
end

function Environment(R::Real, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0))
	return Environment(0..R, ocn, bty, ati)
end

## End of Module
end