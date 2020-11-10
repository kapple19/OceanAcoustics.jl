export Boundary
export Celerity
export Medium
export Environment
export Position
export Signal
export Fan
export Source
export Scenario
export Ray
export Trace
export Beam
export Field
export Grid

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

"""
`Boundary(z, c, ρ) <: OceanAcoustics`

Mutable struct that processes a boundary layer of a vertical slice of an ocean environment.

| Property | Parameter | Type |
|---|---|---|
| `z` | Depth | `Function` |
| `c` | Sound speed | `Function` |
| `ρ` | Density | `Function` |

Each input is versatile. The type passed is used to define a function via linear interpolation. The valid types are as follows.

| Type | Processing |
|---|---|
| `Function` | A univariate function dependent on range. |
| `Tuple{Vr,Vz}`* | A pair of real-valued vectors. The first defines range points, and the second defines respective parameter value. |
| `Real` | A constant boundary parameter value. |

*The full `Tuple{Vr,Vz}` signature is
```
Tuple{Vr, Vz} where {
	Rr <: Real, Vr <: AbstractVector{Rr},
	Rz <: Real, Vz <: AbstractVector{Rz}
}
```
"""
mutable struct Boundary <: OceanAcoustic
	"""
	`z::Function`
	
	Range-dependent boundary depth univariate function.
	"""
	z::Function

	"""
	`callback::ContinuousCallback`
	
	Continuous callback for differential equation solver.
	"""
	callback::ContinuousCallback

	"""
	`Ωr::Interval`
	
	Interval of range values.
	"""
	Ωr::Interval

	"""
	`Ωz::Interval`
	
	Interval of depth values.
	"""
	Ωz::Interval
	
	function Boundary(z::Function)
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
		callback = ContinuousCallback(condition, affect!)
		return new(z, callback)
	end
end

function Boundary(r::Vector{T}, z::Vector{T}) where T <: Real
	zFcn = interpolated_function(r, z)
	bnd = Boundary(zFcn)
	bnd.Ωz = minimum(z)..maximum(z)
	return bnd
end

function Boundary(rz::AbstractArray{T}) where T <: Real
	r = [rng for rng ∈ rz[:, 1]]
	z = [dpt for dpt ∈ rz[:, 2]]
	return Boundary(r, z)
end

function Boundary(z::Real)
	zFcn(r) = z
	bnd = Boundary(zFcn)
	bnd.Ωz = z..z
	return bnd
end

"""
`Celerity <: OceanAcoustic`

struct that processes first and second order partial derivatives of a bivariate celerity function.
"""
struct Celerity <: OceanAcoustic
	c::Function
	∂c_∂r::Function
	∂c_∂z::Function
	∂²c_∂r²::Function
	∂²c_∂r∂z::Function
	∂²c_∂z²::Function
	function Celerity(c::Function)
		c_(x) = c(x[1], x[2])
		∇c_(x) = gradient(c_, x)
		∇c(r, z) = ∇c_([r, z])
		∂c_∂r(r, z) = ∇c(r, z)[1]
		∂c_∂z(r, z) = ∇c(r, z)[2]
	
		∂c_∂r_(x) = ∂c_∂r(x[1], x[2])
		∇∂c_∂r_(x) = gradient(∂c_∂r_, x)
		∇∂c_∂r(r, z) = ∇∂c_∂r_([r, z])
	
		∂c_∂z_(x) = ∂c_∂z(x[1], x[2])
		∇∂c_∂z_(x) = gradient(∂c_∂r_, x)
		∇∂c_∂z(r, z) = ∇∂c_∂z_([r, z])
	
		∂²c_∂r²(r, z) = ∇∂c_∂r(r, z)[1]
		∂²c_∂r∂z(r, z) = ∇∂c_∂r(r, z)[2]
		∂²c_∂z²(r, z) = ∇∂c_∂z(r, z)[2]

		return new(c, ∂c_∂r, ∂c_∂z, ∂²c_∂r², ∂²c_∂r∂z, ∂²c_∂z²)
	end
end

"""
`Medium <: OceanAcoustic`

Mutable struct that processes an acoustic medium in a vertical slice of an ocean environment.

| Property | Parameter | Type |
|---|---|---|
| `c` | Sound Speed | Function |
| `ρ` | Density | Function |
| `α` | Attenuation | Function |

Each input is versatile. The type passed is used to interpolate a function. The valid types are as follows.

| Type | Processing |
|---|---|
| `Function` | A bivariate function dependent on range then depth. |
| `Tuple{Vr,Vz,Mp}`* | Two real-valued vectors followed by a matrix. The first defines range points, the second depth points, then a respective matrix of parameter values. |
| `Tuple{Vr,Tz,Tp}`** | A real-valued vector of ranges followed by two tuples. |
| `Tuple{Vz,Vp}`*** | A pair of real-valued vectors giving depth-parameter value pairs. |
| `Real` | A constant parameter value. |

*The full `Tuple{Vr, Vz, Mp}` signature is
```
Tuple{Vr, Vz, Mp} where {
	Rr <: Real,
	Rz <: Real,
	Rp <: Real,
	Vr <: AbstractVector{Rr},
	Vz <: AbstractVector{Rz},
	Mp <: AbstractArray{Rp}
}
```

**The full `Tuple{Vr, Tz, Tp}` signature is
```
Tuple{Tr, Tz, Tp} where {
	where Nr <: Integer,
	Rr <: Real, Tr <: NTuple{Nr, Rr},
	Rz <: Real, Tz <: NTuple{Nr, AbstractVector{Rz}},
	Rp <: Real, Tp <: NTuple{Nr, AbstractVector{Rp}}
}
```

***The full `Tuple{Vz, Vp}` signature is
```
Tuple{Vz, Vp} where {
	Rz <: Real, Vz <: AbstractVector{Rz},
	Rp <: Real, Vp <: AbstractVector{Rp}
}
```
"""
mutable struct Medium <: OceanAcoustic
	SSP::Celerity
	Ωr::Interval
	Ωz::Interval
	input_data_type::Symbol

	function Medium(c::Function)
		SSP = Celerity(c)
		med = new(SSP)
		med.input_data_type = :Function
		return med
	end
end

function Medium(c::AbstractArray)
	r_ = [rc for rc ∈ c[1, 2:end]]
	z_ = [zc for zc ∈ c[2:end, 1]]
	c_ = c[2:end, 2:end]'
	
	cFcn = interpolated_function(r_, z_, c_)
	med = Medium(cFcn)
	med.input_data_type = :RangeAndDepthDependent
	return med
end

function Medium(z::AbstractVector, c::AbstractVector)
	cMat = vcat([0 0 1.0], hcat(z, c, c))
	med = Medium(cMat)
	med.input_data_type = :DepthDependent
	return med
end

function Medium(c::Real)
	cFcn(r, z) = c
	med = Medium(cFcn)
	med.input_data_type = :Constant
	return med
end

struct Environment <: OceanAcoustic
	Ωr::Interval
	Ωz::Interval
	callback::ContinuousCallback
	ocn::Medium
	bty::Boundary
	ati::Boundary

	function Environment(Ωr::Interval, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0))

		ocn.Ωr = bty.Ωr = ati.Ωr = Ωr

		if !isdefined(bty, :Ωz)
			bty.Ωz = bty.z(Ωr)
		end
		if !isdefined(ati, :Ωz)
			ati.Ωz = ati.z(Ωr)
		end
		ocn.Ωz = Ωz = bty.Ωz ∪ ati.Ωz

		# function c(r, z)
		# 	if ati.z(r) ≤ z ≤ bty.z(r)
		# 		return ocn.SSP.c(r, z)
		# 	else
		# 		return NaN
		# 	end
		# end
		# ocn.SSP = Celerity(c)

		condition(u, t, ray) = (Ωr.hi - Ωr.lo)/2 - abs(u[1] - (Ωr.hi - Ωr.lo)/2)
		affect!(ray) = terminate!(ray)
		callback = ContinuousCallback(condition, affect!)

		return new(Ωr, Ωz, callback, ocn, bty, ati)
	end
end

Environment(R::Real, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0)) = Environment(0..R, ocn, bty, ati)

struct Position <: OceanAcoustic
	r::Real
	z::Real
end

struct Signal <: OceanAcoustic
	f::Real
	ω::Real
	function Signal(f::Real)
		if f ≤ 0
			DomainError(f, "Frequency must be positive.") |> throw
		end
		return new(f, 2π*f)
	end
end

struct Fan
	θ₀s::AbstractVector{R} where R <: Real
	δθ₀s::AbstractVector{R} where R <: Real

	function Fan(θ₀s::AbstractVector{R}) where R <: Real
		if ~all(-π/2 .≤ θ₀s .≤ π/2)
			error("Invalid fan angles.")
		end
		δθ₀s = diff(θ₀s)
		push!(δθ₀s, δθ₀s[end])
		return new(θ₀s, δθ₀s)
	end
end

Fan(θ₀::Real) = Fan([θ₀])

struct Source <: OceanAcoustic
	pos::Position
	sig::Signal
	fan::Fan
end

struct Scenario <: OceanAcoustic
	env::Environment
	src::Source
	name::AbstractString
end

function ray_propagation(scn::Scenario)
	c(r, z) = scn.env.ocn.SSP.c(r, z)
	∂c_∂r(r, z) = scn.env.ocn.SSP.∂c_∂r(r, z)
	∂c_∂z(r, z) = scn.env.ocn.SSP.∂c_∂z(r, z)
	∂²c_∂r²(r, z) = scn.env.ocn.SSP.∂²c_∂r²(r, z)
	∂²c_∂z²(r, z) = scn.env.ocn.SSP.∂²c_∂z²(r, z)
	∂²c_∂r∂z(r, z) = scn.env.ocn.SSP.∂²c_∂r∂z(r, z)

	function ray_propagation!(du, u, p, s)
		r = u[1]
		z = u[2]
		ξ = u[3]
		ζ = u[4]
		τ = u[5]
		pʳ = u[6]
		pⁱ = u[7]
		qʳ = u[8]
		qⁱ = u[9]

		∂²c_∂n²(r, z) = c(r, z)^2 * (
			∂²c_∂r²(r, z) * ζ^2
			- 2∂²c_∂r∂z(r, z) * ξ * ζ
			+ ∂²c_∂z²(r, z) * ξ^2
		)

		du[1] = dr_ds = c(r, z) * ξ
		du[2] = dz_ds = c(r, z) * ζ
		du[3] = dξ_ds = -∂c_∂r(r, z) / c(r, z)^2
		du[4] = dζ_ds = -∂c_∂z(r, z) / c(r, z)^2
		du[5] = dτ_ds = 1/c(r, z)
		du[6] = dpʳ_ds = ∂²c_∂n²(r, z) / c(r, z)^2 * qʳ
		du[7] = dpⁱ_ds = ∂²c_∂n²(r, z) / c(r, z)^2 * qⁱ
		du[8] = dqʳ_ds = c(r, z) * pʳ
		du[9] = dqⁱ_ds = c(r, z) * pⁱ
	end

	callbacks = CallbackSet(
		scn.env.callback,
		scn.env.bty.callback,
		scn.env.ati.callback
	)

	δθ₀s = Vector{Real}(undef, 0)

	r₀ = scn.src.pos.r
	z₀ = scn.src.pos.z
	τ₀ = 0.0

	λ₀ = c(r₀, z₀) / scn.src.sig.f
	pʳ₀ = 1.0
	pⁱ₀ = 0.0
	W₀ = 100λ₀
	qʳ₀ = 0.0
	qⁱ₀ = scn.src.sig.ω * W₀^2 / 2

	TLmax = 100.0
	S = 10^(TLmax/10)
	sSpan = (0.0, S)

	sols = Vector{ODECompositeSolution}(undef, 0)
	for (nRay, θ₀) ∈ enumerate(scn.src.fan.θ₀s)
		push!(δθ₀s, scn.src.fan.δθ₀s[nRay])

		ξ₀ = cos(θ₀) / c(r₀, z₀)
		ζ₀ = sin(θ₀) / c(r₀, z₀)
		u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, pʳ₀, pⁱ₀, qʳ₀, qⁱ₀]
	
		prob = ODEProblem(ray_propagation!, u₀, sSpan)
		push!(
			sols,
			solve(
				prob,
				AutoVern7(Rodas4()),
				callback = callbacks,
				reltol=1e-8, abstol=1e-8
			)
		)
	end

	return sols, δθ₀s
end

struct Ray <: OceanAcoustic
	Ωs::Interval
	δθ₀::Real
	r::Function
	z::Function
	ξ::Function
	ζ::Function
	τ::Function
	p::Function
	q::Function
	θ::Function
	c::Function

	function Ray(sol::ODECompositeSolution, δθ₀::Real)
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
		return new(0..S, δθ₀, r, z, ξ, ζ, τ, p, q, θ, c)
	end
end

struct Trace <: OceanAcoustic
	scn::Scenario
	rays::AbstractVector{R} where R <: Ray

	function Trace(scn::Scenario)
		sols, δθ₀s = ray_propagation(scn)
		rays = Ray.(sols, δθ₀s)
		return new(scn, rays)
	end
end

struct Beam <: OceanAcoustic
	ray::Ray
	p::Function
	function Beam(ray::Ray, src::Source)
		r(s) = ray.r(s)
		z(s) = ray.z(s)
		τ(s) = ray.τ(s)
		p(s) = ray.p(s)
		q(s) = ray.q(s)
		c(s) = ray.c(s)

		c₀ = c(0)
		ω = src.sig.ω
		λ₀ = c₀/ω
		q₀ = q(0)
		θ₀ = ray.θ(0)
		δθ₀ = ray.δθ₀
		
		A = δθ₀/c₀ * exp(im*π/4) * √(q₀ * ω * cos(θ₀) / 2π)
		pressure_beam(s, n) = A * √(c(s) / r(s) / q(s)) * exp(-im * ω * (τ(s) + p(s)/q(s) * n^2 / 2))

		return new(ray, pressure_beam)
	end
end

struct Coherence <: OceanAcoustic
	coherence::Symbol
	pre_summation::Function
	post_summation::Function
end

function Coherence(coh::Symbol = :Coherent)
	if coh == :Coherent
		return Coherence(
			coh,
			p -> p,
			p -> p
		)
	elseif coh == :Incoherent
		return Coherence(
			coh,
			p -> abs(p),
			p -> p
		)
	else
		error("Unrecognized coherence.")
	end
end

Coherence(coh::AbstractString) = Symbol(coh) |> Coherence

struct Field <: OceanAcoustic
	scn::Scenario
	beams::AbstractVector{Beam}
	coh::Coherence
	p::Function

	function Field(
		scn::Scenario,
		beams::AbstractVector{B},
		coh::Coherence = Coherence()
		) where B <: Beam

		function pressure(r, z)
			p = [
				beam.p(s, n)
				for beam ∈ beams
				for (s, n) ∈ closest_points(r, z, beam)
			] |> sum
			# p = complex(0.0)
			# for beam ∈ beams, (s, n) ∈ closest_points(r, z, beam)
			# 	p += beam.p(s, n)
			# end
		end

		return new(
			scn,
			beams,
			coh,
			pressure
		)
	end
end

Field(trc::Trace, coh::Coherence = Coherence()) = Field(trc.scn, Beam.(trc.rays, trc.scn.src), coh)

Field(scn::Scenario, coh::Coherence = Coherence()) = Field(Trace(scn), coh)

struct Sonar <: OceanAcoustic
	
end

using Base.Threads

struct Grid <: OceanAcoustic
	scn::Scenario
	r::AbstractVector{Rr} where Rr <: Real
	z::AbstractVector{Rz} where Rz <: Real
	p::AbstractArray{Cp, 2} where Cp <: Number
	TL::AbstractArray{RTL, 2} where RTL <: Real

	function Grid(
		fld::Field,
		r::AbstractVector{Rr},
		z::AbstractVector{Rz}
		) where {Rr <: Real, Rz <: Real}
		
		DEF_NAME = "Pressure Grid"
		progress_name(name) = length(name) == 0 ? DEF_NAME : name * ": " * DEF_NAME * " "
		pn = progress_name(fld.scn.name)

		function pressure(r, z)
			if fld.scn.env.ati.z(r) < z < fld.scn.env.bty.z(r)
				return fld.p(r, z)
			else
				return NaN
			end
		end

		# p .= @showprogress 1 pn [pressure(r′, z′) for z′ ∈ z, r′ ∈ r]

		Nr = length(r)
		Nz = length(z)
		rGrid = reshape(
			r' .* ones(Nz, Nr),
			(Nr*Nz,)
		)
		zGrid = reshape(
			z .* ones(Nz, Nr),
			(Nr*Nz,)
		)

		p′ = zeros(Complex, Nr*Nz)
		@show nthreads()
		@time @threads for n ∈ 1:Nr*Nz
		# @time @showprogress 1 pn for n ∈ 1:Nr*Nz
			# r′ = rGrid[n]
			# z′ = zGrid[n]
			@inbounds p′[n] = pressure(rGrid[n], zGrid[n])
		end
		p = reshape(p′, Nz, Nr)

		TL = min.(100, SonarEqs.transmission_loss.(p))

		return new(fld.scn, r, z, p, TL)
	end
end

function Grid(fld::Field, Nr::Integer, Nz::Integer)
	r = gridpoints(fld.scn.env.Ωr, Nr)
	z = gridpoints(fld.scn.env.Ωz, Nz)

	return Grid(fld, r, z)
end

function Grid(scn::Scenario, args...)
	fld = Field(scn)
	return Grid(fld, args...)
end

Grid(oac::OceanAcoustic) = Grid(oac, 51, 45)
