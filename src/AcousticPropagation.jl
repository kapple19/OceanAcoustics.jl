# Dependencies
using OrdinaryDiffEq:
ContinuousCallback,
CallbackSet,
ODEProblem,
solve,
terminate!,
ODESolution,
Rodas4,
AutoVern7
using ForwardDiff: derivative
using Base: broadcastable
using Roots: find_zeros
using IntervalArithmetic:
(..),
Interval

# Exports:
# * Environment
export Boundary
export Medium
export Environment

# * Scenario
export Position
export Signal
export Source
export Scenario

# * Propagation
export Ray
export Beam
export Field

"""
	boundary_reflection(t_inc::Vector, t_bnd::Vector) -> t_rfl::Vector

Calculates the reflection ray tangent vector `r_rfl` for an incident ray tangent vector `t_inc` reflecting against a boundary with tangent vector `t_bnd`.
"""
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
	dz_dr::Function
	condition::Function
	affect!::Function
	R::Real
	
	"""
		Boundary(z::Function)

	An ocean boundary storing its depth `z` (metres) as univariate function of range (metres).
	"""
	function Boundary(z::Function, R::Real)
		dz_dr(r) = derivative(z, r)
		condition(u, t, ray) = z(u[1]) - u[2]
		function affect!(ray)
			ξ, ζ = boundary_reflection([ray.u[3], ray.u[4]], [1, dz_dr(ray.u[1])])
			if ξ < 0
				return terminate!(ray)
			else
				function reflect!(ray)
					ray.u[3] = ξ
					ray.u[4] = ζ
				end
				return reflect!(ray)
			end
		end
		return new(z, dz_dr, condition, affect!, R)
	end
end

"""
	Boundary(r::Vector, z::Vector)

An ocean boundary storing its depth `z` (metres) at range `r` (metres).

The inputted values are interpolated into and stored as a function.
"""
function Boundary(r::Vector{T}, z::Vector{T}) where T <: Real
	zFcn = interpolated_function(r, z)
	return Boundary(zFcn)
end

"""
	Boundary(rz::AbstractArray)

An ocean boundary storing its depth and range as a two-column matrix. The first column contains range (metres), the second column contains the respective depth (metres).

The inputted values are interpolated into and stored as a function.
"""
function Boundary(rz::AbstractArray{T}) where T <: Real
	r = [rng for rng ∈ rz[:, 1]]
	z = [dpt for dpt ∈ rz[:, 2]]
	return Boundary(r, z)
end

"""
	Boundary(z::Real)

An ocean boundary storing its depth `z` (metres) as a constant.

The inputted valued is interpolated into and stored as a function.
"""
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
	∂²c_∂z²::Function
end

function Celerity(c::Function)
	∂c_∂r, ∂c_∂z, ∂²c_∂r², ∂²c_∂r∂z, ∂²c_∂z² = bivariate_partial_derivatives(c)
	return Celerity(c, ∂c_∂r, ∂c_∂z, ∂²c_∂r², ∂²c_∂r∂z, ∂²c_∂z²)
end

function Celerity(z::AbstractVector, c::AbstractVector)
	cFcn = interpolated_function(z, c)
	return Celerity(cMat)
end

function Celerity(c::Real)
	cFcn(r, z) = c
	return Celerity(cFcn)
end

struct Medium <: OceanAcoustic
	SSPₚ::Celerity
	SSPₛ::Celerity
	# ρ
	# α
	# T
end

function Medium(c::Union{Function, Real})
	SSPₚ = Celerity(c)
	return Medium(SSPₚ)
end

function fluid_medium(SSPₚ::Celerity)
	SSPₛ = Celerity((r, z) -> 0.0)
	return Medium(SSPₚ, SSPₛ)
end

function solid_medium(SSPₚ::Celerity, SSPₛ::Celerity)
	return Medium(SSPₚ, SSPₛ)
end

struct Environment <: OceanAcoustic
	media::AbstractVector{M} where M <: Medium
	boundaries::AbstractVector{B} where B <: Boundary
	Ωrange::Interval
	Ωdepth::Interval
	Ωdepths::Union{RANGE, AbstractVector{RANGE}} where RANGE <: Interval
	c

	function Environment(
		media::AbstractVector{M},
		boundaries::AbstractVector{B},
		Ωrange::Interval
		) where {M <: Medium, B <: Boundary}

		if length(media) ≠ length(boundaries) + 1
			DimensionMismatch("Environment requires the media on each side of each boundary.") |> throw
		end

		Ωdepths = Vector{Interval}(undef, 0)
		for nBnd = 1:length(boundaries) - 1
			zₙ(r) = boundaries[nBnd].z(r)
			zₙ₊₁(r) = boundaries[nBnd + 1].z(r)
			Δz(r) = zₙ(r) - zₙ₊₁(r)
			Ω_Δz = Δz(Ωrange)
			if Ω_Δz.lo ≤ 0
				nBnd₊ = nBnd + 1
				ErrorException("Boundaries $nBnd and $nBnd₊ should never meet.") |> throw
			end

			Ωl = zₙ(Ωrange)
			Ωu = zₙ₊₁(Ωrange)

			push!(Ωdepths, Ωl.lo..Ωu.hi)

			for celerity ∈ (:SSPₚ, :SSPₛ)
				SSP = getproperty(media[nBnd], celerity)
				cₙ = SSP.c
				function c(r′)
					Ω_z(r) = zₙ(r)..zₙ₊₁(r)
					cₙ(r′, Ω_z(r′))
				end
				Ωc = c(Ωrange)
				if Ωc.lo ≤ 0
					ErrorException("Ocean SSP $nBnd must be positive between its respective boundaries.")
				end
			end
		end
		Ωdepth = Ωdepths[1] ∪ Ωdepths[end]

		function c(r, z)
			for nBnd ∈ 1:length(boundaries) - 1
				zₙ(r) = boundaries[nBnd].z(r)
				zₙ₊₁(r) = boundaries[nBnd + 1].z(r)
				if z == zₙ(r)
					return media[nBnd].SSPₚ.c(r, z)
				elseif z == zₙ₊₁(r)
					return media[nBnd+2].SSPₚ.c(r, z)
				elseif zₙ(r) < z < zₙ₊₁(r)
					return media[nBnd+1].SSPₚ.c(r, z)
				end
			end
		end

		return new(media, boundaries, Ωrange, Ωdepth, Ωdepths, c)
	end
end

function Environment(
	media::AbstractVector{M},
	boundaries::AbstractVector{B},
	R::Real
	) where {M <: Medium, B <: Boundary}
	
	Ωrange = 0..R
	return Environment(media, boundaries, Ωrange)
end

"""
	Position(r::Real, z::Real)

Position in 2D slice of ocean, with range `r` (metres) and depth `z` (metres).
"""
struct Position <: OceanAcoustic
	r::Real
	z::Real
end

"""
	Signal(f::Real)

Parameters for a signal with frequency `f` (Hertz).
"""
struct Signal <: OceanAcoustic
	f::Real

	function Signal(f::Real)
		if f ≤ 0
			DomainError(f, "Frequency must be positive.") |> throw
		end
		return new(f)
	end
end

struct Source <: OceanAcoustic
	θ₀::Union{THETA, AbstractVector{THETA}} where THETA <: Real
	δθ₀::Union{THETA, AbstractVector{THETA}} where THETA <: Real
	pos::Position
	sig::Signal
end

function Source(θ₀::AbstractVector{THETA}, pos::Position, sig::Signal) where THETA <: Real
	δθ₀ = diff(θ₀)
	push!(δθ₀, δθ₀[end])
	return Source(θ₀, δθ₀, pos, sig)
end

function Source(θ₀::Real, pos::Position, sig::Signal)
	δθ₀ = 1.0
	return Source(θ₀, δθ₀, pos, sig)
end

struct Scenario <: OceanAcoustic
	sources::Union{SOURCE, AbstractVector{SOURCE}} where SOURCE <: Source
	env::Environment
	layers::Union{LOCATION, AbstractVector{LOCATION}} where LOCATION <: Integer
end

function propagation_problem(scenario::Scenario)
	# Initialize vector of `ODEProblem`s and `CallbackSet`s
	prop_probs = Vector{ODEProblem}(undef, 0)
	callbacks = Vector{CallbackSet}(undef, 0)
	δθ₀s = Vector{Real}(undef, 0)

	# Range condition
	Ωr = scenario.env.Ωrange
	rng_condition(u, t, ray) = (Ωr.hi - Ωr.lo)/2 - abs(u[1] - (Ωr.hi + Ω.lo)/2)
	rng_affect!(ray) = terminate!(ray)

	# Universal Initial Conditions
	τ₀ = 0.0
	p₀ʳ = 1.0
	p₀ⁱ = 0.0
	q₀ʳ = 0.0

	# Maximum Path Length
	TLmax = 100.0
	S = 10^(TLmax/10.0)
	sSpan = (0., S)

	# Loop on each sound source
	for (nSrc, src) ∈ enumerate(scenario.sources)
		# Get Environment
		nLayer = scenario.layers[nSrc]
		ocn = scenario.env.media[nLayer]
		ati = scenario.env.boundaries[nLayer]
		bty = scenario.env.media[nLayer + 1]
		Ωz = scenario.env.Ωdepths[nLayer]

		# Get Celerity
		c = ocn.SSPₚ.c
		∂c_∂r = ocn.SSPₚ.∂c_∂r
		∂c_∂z = ocn.SSPₚ.∂c_∂z
		
		# Equations:
		# * eikonal equation
		# * transport equation (first order)
		# * dynamic ray equations
		function propagation!(du, u, p, s)
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
				ocn.SSPₚ.∂²c_∂r²(r, z)*ζ^2
				- 2ocn.SSPₚ.∂²c_∂r∂z(r, z)*ξ*ζ
				+ ocn.SSPₚ.∂²c_∂z²(r, z)*ξ^2
			)

			# Differential Equations:
			# * Eikonal
			du[1] = dr_ds = c(r, z)*ξ
			du[2] = dz_ds = c(r, z)*ζ
			du[3] = dξ_ds = -∂c_∂r(r, z)/c(r, z)^2
			du[4] = dζ_ds = -∂c_∂z(r, z)/c(r, z)^2

			# * Time
			du[5] = dτ_ds = 1/c(r, z)

			# * Dynamic Ray Equations
			du[6] = dpʳ_ds = ∂²c_∂n²(r, z)/c(r, z)^2*qʳ
			du[7] = dpⁱ_ds = ∂²c_∂n²(r, z)/c(r, z)^2*qⁱ
			du[8] = dqʳ_ds = c(r, z)*pʳ
			du[9] = dqⁱ_ds = c(r, z)*pⁱ
		end

		# Boundary Callbacks
		CbRng = ContinuousCallback(rng_condition, rng_affect!)
		CbBty = ContinuousCallback(bty.condition, bty.affect!)
		CbAti = ContinuousCallback(ati.condition, ati.affect!)

		# Callback Set
		push!(callbacks, CallbackSet(CbRng, CbBty, CbAti))

		# Initial conditions for source
		r₀ = src.pos.r
		z₀ = src.pos.z
		ω = src.sig.f
		λ₀ = c(r₀, z₀)/src.sig.f
		W₀ = 100λ₀ # 10..50
		q₀ⁱ = ω*W₀^2/2

		# Loop for each angle
		for (θ₀, δθ₀) ∈ [(src.θ₀[nθ₀], src.δθ₀[nθ₀]) for nθ₀ ∈ eachindex(src.θ₀)]
			# Initial condition for angle
			ξ₀ = cos(θ₀)/c(r₀, z₀)
			ζ₀ = sin(θ₀)/c(r₀, z₀)
			push!(δθ₀s, src.δθ₀)

			# Initial Condition
			u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀ʳ, p₀ⁱ, q₀ʳ, q₀ⁱ]

			# Accumulate ODEProblems
			push!(prop_probs, ODEProblem(propagation!, u₀, sSpan))
		end
	end
	return prop_probs, callbacks, δθ₀s
end

function solve_propagation(prop_prob::ODEProblem, callback::CallbackSet)
	RaySol = solve(prop_prob, AutoVern7(Rodas4()), callback = callback, reltol=1e-8, abstol=1e-8)
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
end

function Ray(RaySol::ODESolution, δθ₀::Real)
	S = RaySol.t[end]
	r(s) = RaySol(s, idxs = 1)
	z(s) = RaySol(s, idxs = 2)
	ξ(s) = RaySol(s, idxs = 3)
	ζ(s) = RaySol(s, idxs = 4)
	τ(s) = RaySol(s, idxs = 5)
	p(s) = RaySol(s, idxs = 6) + im*RaySol(s, idxs = 7)
	q(s) = RaySol(s, idxs = 8) + im*RaySol(s, idxs = 9)
	θ(s) = atan(ζ(s)/ξ(s))
	c(s) = cos(θ(s))/ξ(s)

	return Ray(S, r, z, ξ, ζ, τ, p, q, θ, c, δθ₀)
end

function Ray(scenario::Scenario)
	prop_probs, callbacks, δθ₀s = propagation_problem(scenario)
	
	RaySols = solve_propagation.(Prob, CbBnd)
	
	return Ray.(RaySols, δθ₀s)
end

struct Beam <: OceanAcoustic
	ray
	b::Function
	W::Function
end

function Beam(ray::Ray)
	r(s) = ray.r(s)
	z(s) = ray.z(s)
	τ(s) = ray.τ(s)
	p(s) = ray.p(s)
	q(s) = ray.q(s)
	c(s) = ray.c(s)
	W(s) = sqrt(-2/ω/imag(p(s)/q(s)))

	c₀ = c(0)
	ω = 2π*src.sig.f
	λ₀ = c₀/src.sig.f
	W₀ = W(0)
	q₀ = q(0)
	δθ₀ = ray.δθ₀

	A = δθ₀/c₀ * exp(im*π/4)*sqrt(q₀*ω*cos(θ₀)/2π)
	b(s, n) = A * sqrt(c(s)/r(s)/q(s)) * exp(-im*ω * (τ(s) + p(s)/q(s)*n^2/2))

	return Beam(ray, b, W)
end

function Beam(scenario::Scenario)
	rays = Ray(scenario)
	return Beam.(rays)
end

function closest_points(r, z, beam)
	Q(s) = (beam.ray.r(s) - r)^2 + (beam.ray.z(s) - z)^2
	dQ(s) = derivative(Q, s)
	sMins = find_zeros(dQ, 0, beam.ray.S)
	d²Q(s) = derivative(dQ, s)
	# min_cond(s) = d²Q(s) > 0 && beam.W(s) > sqrt(Q(s))
	min_cond(s) = d²Q(s) > 0
	min_cond.(sMins)
	filter!(min_cond, sMins)
	return sMins, sqrt.(Q.(sMins))
end

function add_to_pressure(r::Real, z::Real, beam::Beam, δθ₀::Real, coh_pre::Function)
	sMins, nMins = closest_points(r, z, beam)
	p = complex(0)
	for (n, sMin) ∈ enumerate(sMins)
		p += coh_pre(δθ₀ * beam.b(sMin, nMins[n]))
	end
	return p
end

struct Field <: OceanAcoustic
	beams::AbstractVector{Beam}
	src::Source
	ocn::Medium
	bty::Boundary
	ati::Boundary
	p::Function
	TL::Function
end

function Field(beams::AbstractVector{T}, src::Source, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0.0)) where T <: Beam

	coh_pre(p) = p
	coh_post(p) = p

	function pressure(r::Real, z::Real)
		p = complex(0.0)
		for (n, beam) ∈ enumerate(beams)
			p += add_to_pressure(r, z, beam, δθ₀[n], coh_pre)
		end
		return coh_post(p)
		
	end

	function transmission_loss(r::Real, z::Real)
		pAbs = abs(pressure(r, z))
		if isnan(pAbs)
			return 0.0
		else
			TL = max(0.0, min(100.0, -20log10(pAbs)))
			return TL
		end
	end

	return Field(beams, src, ocn, bty, ati, pressure, transmission_loss)
end

function Field(θ₀s::AbstractVector{T}, src::Source, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0.0)) where T <: Real
	beams = Beam.(θ₀s, src, ocn, bty, ati)
	return Field(beams, src, ocn, bty, ati)
end
