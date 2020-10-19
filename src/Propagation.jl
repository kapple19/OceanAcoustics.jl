using ForwardDiff: derivative
using IntervalArithmetic: Interval, (..)
using OrdinaryDiffEq: ContinuousCallback

export (..)
export Boundary
export Medium
export Environment

const SOUND_SPEED_AIR = 343

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
	z::Function
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

function Boundary(z::Any, c::Any)
	zFcn = univariate_interpolation(z...)
	cFcn = univariate_interpolation(c...)
	return Boundary(zFcn, cFcn)
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
		∂c_∂r, ∂c_∂z, ∂²c_∂r², ∂²c_∂r∂z, ∂²c_∂z∂r, ∂²c_∂z² = bivariate_derivatives(c)
		return new(
			c, ∂c_∂r, ∂c_∂z,
			∂²c_∂r²,
			∂²c_∂r∂z, ∂²c_∂z∂r,
			∂²c_∂z²)
	end
end

struct Medium <: OceanAcoustic
	SSPₚ::Celerity
	SSPₛ::Celerity
	
	function Medium(
		SSPₚ::Celerity,
		SSPₛ::Celerity = Celerity(0)
		)
		return new(SSPₚ, SSPₛ)
	end
end

function Medium(cₚ::Function, cₛ::Function)
	SSPₚ = Celerity(cₚ)
	SSPₛ = Celerity(cₛ)
	return Medium(SSPₚ, SSPₛ)
end

function Medium(cₚ::Any, cₛ::Any = 0)
	cₚFcn = bivariate_interpolation(cₚ...)
	cₛFcn = bivariate_interpolation(cₛ...)
	return Medium(cₚFcn, cₛFcn)
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
		ati::Boundary = Boundary(0, SOUND_SPEED_AIR))

		Ωz_ati = ati.z(Ωr)
		Ωz_bty = bty.z(Ωr)
		Ωz = Ωz_ati ∪ Ωz_bty

		return new(Ωr, Ωz, ocn, bty, ati)
	end
end

function Environment(R::Real, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0, SOUND_SPEED_AIR))
	return Environment(0..R, ocn, bty, ati)
end