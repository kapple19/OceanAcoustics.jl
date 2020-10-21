export (..)
export Boundary
export Medium
export Environment
export Position
export Signal
export Spark
export Fan
export Source
export Scenario
export Ray
export Trace

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
				terminate!(ray)
			else
				function reflect!(ray)
					ray.u[3] = ξₒ
					ray.u[4] = ζₒ
				end
				reflect!(ray)
			end
		end
		callback = ContinuousCallback(condition, affect!)

		return new(z, c, callback)
	end
end

function Boundary(z::Any, c::Any)
	zFcn = univariate_interpolation(z)
	cFcn = univariate_interpolation(c)
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
	SSP::Celerity
	Medium(SSP::Celerity) = new(SSP)
end

function Medium(c::Function)
	SSP = Celerity(c)
	return Medium(SSP)
end

function Medium(c::Tuple)
	cFcn = bivariate_interpolation(c...)
	return Medium(cFcn)
end

struct Environment <: OceanAcoustic
	Ωr::Interval
	Ωz::Interval
	callback::ContinuousCallback
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

		function parse_interval(Ω)
			if Ω isa Interval
				return Ω
			elseif length(Ω) == 2
				return Ω = Ω[1]..Ω[2]
			elseif length(Ω) == 1
				return Ω = Ω..Ω
			else
				error("Unknown parsing error.")
			end
		end
		Ωz_ati = parse_interval(Ωz_ati)
		Ωz_bty = parse_interval(Ωz_bty)
		Ωz = Ωz_ati ∪ Ωz_bty
		Ωz = parse_interval(Ωz)

		condition(u, t, ray) = (Ωr.hi - Ωr.lo)/2 - abs(u[1] - (Ωr.hi + Ωr.lo)/2)
		affect!(ray) = terminate!(ray)
		callback = ContinuousCallback(condition, affect!)

		return new(Ωr, Ωz, callback, ocn, bty, ati)
	end
end

function Environment(R::Real, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0, SOUND_SPEED_AIR))
	return Environment(0..R, ocn, bty, ati)
end

struct Position <: OceanAcoustic
	r::Real
	z::Real
	Position(r::Real, z::Real) = new(r, z)
end

struct Signal <: OceanAcoustic
	f::Real
	ω::Real
	Signal(f::Real) = f ≤ 0 ? error("Signal frequency must be positive.") : return new(f, 2π*f)
end

struct Spark <: OceanAcoustic
	θ₀::Real
	Spark(θ₀::Real) = new(θ₀)
end

struct Fan <: OceanAcoustic
	spks::AbstractVector{S} where S <: Spark
	δθ₀s::AbstractVector{R} where R <: Real
	function Fan(
		spks::AbstractVector{S},
		δθ₀::AbstractVector{R}
		) where {S <: Spark, R <: Real}
		return new(spks, δθ₀)
	end
end

function Fan(spks::AbstractVector{S}) where S <: Spark
	δθ₀s = [spks[nSpk].θ₀ for nSpk ∈ eachindex(spks)] |> diff
	push!(δθ₀s, δθ₀s[end])
	return Fan(spks, δθ₀s)
end

function Fan(spk::Spark)
	δθ₀ = 1.0
	return Fan([spk], [δθ₀])
end

function Fan(θ₀s::AbstractVector{R}) where R <: Real
	return Fan(Spark.(θ₀s))
end

function Fan(θ₀::Real)
	return Fan(Spark(θ₀))
end

struct Source <: OceanAcoustic
	pos::Position
	sig::Signal
	fan::Fan
end

struct Scenario <: OceanAcoustic
	env::Environment
	srcs::Union{S, AbstractVector{S}} where S <: Source
	name::AbstractString
	function Scenario(
		env::Environment,
		srcs::AbstractVector{S},
		name::AbstractString = "") where S <: Source
		return new(env, srcs, name)
	end
end

Scenario(env::Environment, src::Source, name::AbstractString = "") = Scenario(env, [src], name)

function propagation(sno::Scenario)
	callbacks = CallbackSet(
		sno.env.callback,
		sno.env.ati.callback,
		sno.env.bty.callback
	)

	c(r, z) = sno.env.ocn.SSP.c(r, z)
	∂c_∂r(r, z) = sno.env.ocn.SSP.∂c_∂r(r, z)
	∂c_∂z(r, z) = sno.env.ocn.SSP.∂c_∂z(r, z)
	∂²c_∂r²(r, z) = sno.env.ocn.SSP.∂²c_∂r²(r, z)
	∂²c_∂r∂z(r, z) = sno.env.ocn.SSP.∂²c_∂r∂z(r, z)
	∂²c_∂z²(r, z) = sno.env.ocn.SSP.∂²c_∂z²(r, z)
	
	function propagate!(du, u, p, s)
		r = u[1]
		z = u[2]
		ξ = u[3]
		ζ = u[4]
		τ = u[5]
		pʳ = u[6]
		pⁱ = u[7]
		qʳ = u[8]
		qⁱ = u[9]

		# Second partial derivative WRT ray normal
		∂²c_∂n²(r, z) = c(r, z)^2*(
			∂²c_∂r²(r, z)*ζ^2
			- 2∂²c_∂r∂z(r, z)*ξ*ζ
			+ ∂²c_∂z²(r, z)*ξ^2
		)

		# Differential Equations:
		# * Eikonal Equation (ray-parameterized)
		du[1] = dr_ds = c(r, z)*ξ
		du[2] = dz_ds = c(r, z)*ζ
		du[3] = dξ_ds = -∂c_∂r(r, z)/c(r, z)^2
		du[4] = dζ_ds = -∂c_∂z(r, z)/c(r, z)^2

		# * Transport (first order)
		du[5] = dτ_ds = 1/c(r, z)

		# * Dynamic Ray Equations (split into real and imaginary parts)
		du[6] = dpʳ_ds = ∂²c_∂n²(r, z)/c(r, z)^2*qʳ
		du[7] = dpⁱ_ds = ∂²c_∂n²(r, z)/c(r, z)^2*qⁱ
		du[8] = dqʳ_ds = c(r, z)*pʳ
		du[9] = dqⁱ_ds = c(r, z)*pⁱ
	end

	# Persistent initial conditions
	τ₀ = 0.0
	p₀ʳ = 1.0
	p₀ⁱ = 0.0
	q₀ʳ = 0.0

	# Maximum ray length
	TLmax = 100.0
	S = 10^(TLmax/10)
	sspan = (0.0, S)

	# Initialize vector of ODEProblems and ray launch deltas
	probs = Vector{ODEProblem}(undef, 0)
	δθ₀s = Vector{Real}(undef, 0)
	for src ∈ sno.srcs
		# Initial conditions for each source
		r₀ = src.pos.r
		z₀ = src.pos.z
		λ₀ = c(r₀, z₀)/src.sig.f
		W₀ = 100λ₀ # 10..50
		q₀ⁱ = src.sig.ω * W₀^2/2

		for (nSpk, spk) ∈ enumerate(src.fan.spks)
			# Populate launch angle deltas
			push!(δθ₀s, src.fan.δθ₀s[nSpk])

			# Initial conditions for each angle
			ξ₀ = cos(spk.θ₀)/c(r₀, z₀)
			ζ₀ = sin(spk.θ₀)/c(r₀, z₀)

			# Initial Condition
			u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀ʳ, p₀ⁱ, q₀ʳ, q₀ⁱ]

			# Accumulate ODEProblems
			push!(probs, ODEProblem(propagate!, u₀, sspan))
		end
	end

	function propagate(prob::ODEProblem)
		sol = @time solve(
			prob,
			AutoVern7(Rodas4()),
			callback = callbacks
		)
	end

	return propagate.(probs), δθ₀s
end

struct Ray <: OceanAcoustic
	S::Real
	r::Function
	z::Function
	ξ::Function
	ζ::Function
	τ::Function
	p::Function
	q::Function
	θ::Function
	c::Function
	δθ₀::Real
	sol::AbstractODESolution
	function Ray(sol::AbstractODESolution, δθ₀)
		S = sol.t[end]
		r(s) = sol(s, idxs = 1)
		z(s) = sol(s, idxs = 2)
		ξ(s) = sol(s, idxs = 3)
		ζ(s) = sol(s, idxs = 4)
		τ(s) = sol(s, idxs = 5)
		p(s) = sol(s, idxs = 6) + im*sol(s, idxs = 7)
		q(s) = sol(s, idxs = 8) + im*sol(s, idxs = 9)
		θ(s) = atan(ζ(s)/ξ(s))
		c(s) = cos(θ(s))/ξ(s)
		return new(S, r, z, ξ, ζ, τ, p, q, θ, c, δθ₀, sol)
	end
end

struct Trace <: OceanAcoustic
	sno::Scenario
	rays::AbstractVector{R} where R <: Ray
	function Trace(sno::Scenario)
		sols, δθ₀s = propagation(sno)
		return new(sno, Ray.(sols, δθ₀s))
	end
end
