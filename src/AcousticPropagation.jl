using Interpolations:
LinearInterpolation,
Flat
using OrdinaryDiffEq:
ContinuousCallback,
CallbackSet,
ODEProblem,
solve,
terminate!,
ODESolution,
Rodas4,
AutoVern7
using ForwardDiff:
derivative,
gradient
using Base: broadcastable
using Roots: find_zeros
using GRUtils:
Figure,
plot!,
contourf!,
hold!,
yflip!,
xlabel!,
ylabel!,
colorscheme!,
color,
title!,
gcf,
display

export Position
export Signal
export Source
export Boundary
export Medium
export Ray
export Beam
export Receiver
export Field
export acoustic_plot!
export acoustic_plot

function interpolated_function(x, y)
	Itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
	return ItpFcn(x::Real) = Itp(x)
end
function interpolated_function(x, y, z)
	Itp = LinearInterpolation((x, y), z, extrapolation_bc = Flat())
	return ItpFcn(x::Real, y::Real) = Itp(x, y)
end

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

"""
	Position(r::Real, z::Real)

Position in 2D slice of ocean, with range `r` (metres) and depth `z` (metres).
"""
struct Position
	r::Real
	z::Real
end

"""
	Signal(f::Real)

Parameters for a signal with frequency `f` (Hertz).
"""
struct Signal
	f::Real

	function Signal(f::Real)
		if f ≤ 0
			throw(DomainError(f, "Frequency must be positive."))
		end
		return new(f)
	end
end

"""
	Source(pos::Position, sig::Signal)

An ocean sound source with position `pos` and signal `sig`.
"""
struct Source
	pos::Position
	sig::Signal
end

struct Boundary
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
function Boundary(
	r::Vector{T},
	z::Vector{T},
	R::T = r[end]) where T <: Real
	zFcn = interpolated_function(r, z)
	return Boundary(zFcn, R)
end

"""
	Boundary(rz::AbstractArray)

An ocean boundary storing its depth and range as a two-column matrix. The first column contains range (metres), the second column contains the respective depth (metres).

The inputted values are interpolated into and stored as a function.
"""
function Boundary(rz::AbstractArray{T}, R::T = rz[end, 1]) where T <: Real
	r = [rng for rng ∈ rz[:, 1]]
	z = [dpt for dpt ∈ rz[:, 2]]
	return Boundary(r, z, R)
end

"""
	Boundary(z::Real)

An ocean boundary storing its depth `z` (metres) as a constant.

The inputted valued is interpolated into and stored as a function.
"""
function Boundary(z::Real, R::Real = 1.0)
	zFcn(r) = z
	return Boundary(zFcn, R)
end

struct Medium
	c::Function
	∂c_∂r::Function
	∂c_∂z::Function
	∂²c_∂r²::Function
	∂²c_∂r∂z::Function
	∂²c_∂z²::Function
	R::Real
	Z::Real

	"""
		Medium(c::Function, R::Real, Z::Real)
	
	An acoustic medium storing the sound speed `c` (m/s) as a bivariate function of range and depth, with a maximum range `R` (metres) and maximum depth `Z`.

	The following derivatives are also computed and stored as bivariate functions:
	* `∂c_∂r(r, z)`: ∂c/∂r
	* `∂c_∂z(r, z)`: ∂c/∂z
	* `∂²c_∂r²(r, z)`: ∂²c/∂r²
	* `∂²c_∂r∂z(r, z)`: ∂²c/∂r∂z
	* `∂²c_∂z²(r, z)`: ∂²c/∂z²
	"""
	function Medium(c::Function, R::Real, Z::Real)
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
	
		return new(c, ∂c_∂r, ∂c_∂z, ∂²c_∂r², ∂²c_∂r∂z, ∂²c_∂z², R, Z)
	end
end

"""
	Medium(c::AbstractArray, R::Real = c[end, 1], Z::Real = c[1, end])

An acoustic medium storing the sound speed `c` as an array with values `2:end` in the first row as range (metres), values `2:end` in the first column as depth (metres) and values `[2:end, 2:end]` as the respective sound speed (m/s) grid.

The medium maximal range `R` (metres) and maximal depth `Z` are also stored. Their defaults are their respective last values of the given range/depth in the inputted sound speed array.

The inputted values are interpolated into and stored as a bivariate function of range and depth.

The following derivatives are also computed and stored:
* `∂c_∂r(r, z)`: ∂c/∂r
* `∂c_∂z(r, z)`: ∂c/∂z
* `∂²c_∂r²(r, z)`: ∂²c/∂r²
* `∂²c_∂r∂z(r, z)`: ∂²c/∂r∂z
* `∂²c_∂z²(r, z)`: ∂²c/∂z²
"""
function Medium(c::AbstractArray, R::Real = c[end, 1], Z::Real = c[1, end])
	r_ = [rc for rc ∈ c[1, 2:end]]
	z_ = [zc for zc ∈ c[2:end, 1]]
	c_ = c[2:end, 2:end]'
	
	cFcn = interpolated_function(r_, z_, c_)
	return Medium(cFcn, R, Z)
end

"""
	Medium(z::AbstractVector, c::AbstractVector, R::Real, Z::Real = z[end])

An acoustic medium storing the sound speed `c` (m/s) as a vector with corresponding depths `z` (metres).

The medium maximal range `R` (metres) and medium maximal depth `Z` (metres) are also stored. The default of `Z` is the last value in `z`.

The sound speed grid values are interpolated into and stored as a bivariate function of range `r` (metres) and depth `z` (metres) as `c(r, z)`.

The following derivatives are also computed and stored:
* `∂c_∂r(r, z)`: ∂c/∂r
* `∂c_∂z(r, z)`: ∂c/∂z
* `∂²c_∂r²(r, z)`: ∂²c/∂r²
* `∂²c_∂r∂z(r, z)`: ∂²c/∂r∂z
* `∂²c_∂z²(r, z)`: ∂²c/∂z²

Note that for this dispatch, all range derivatives are zero due to range-independence.
"""
function Medium(z::AbstractVector, c::AbstractVector, R::Real, Z::Real = z[end])
	cMat = vcat([0 0 R], hcat(z, c, c))
	return Medium(cMat, R, Z)
end

"""
	Medium

An acoustic medium storing sound speed `c` (m/s) as a constant, along with the maximal medium range `R` and maximal medium depth `Z`.

The sound speed is interpolated and stored as a function. Derivatives are also calculated, but in the case of this dispatch, are zero.
"""
function Medium(c::Real, R::Real, Z::Real)
	cFcn(r, z) = c
	return Medium(cFcn, R, Z)
end

"""
	propagation_problem(θ₀::Real, src::Source, ocn::Medium, bty::Boundary, ati::Boundary) -> prob::ODEProblem, CbBnd::ContinuousCallback

Defines the differential equation problem `Prob` and continuous callback `Cb` for a combination of the eikonal, transport, and dynamic ray equations for a source `src` in an ocean medium `ocn` bounded above by the altimetry `ati` and below by the bathymetry `bty`, for the initial angle `θ₀` (radians) of the ray launched from the specified source position.

The `ODEProblem` `prob` is a type belonging to the `DifferentialEquations.jl` package.

The continuous callback `CbBnd` is defined as condition-affect pairs for checking interactions with the ocean boundaries (in range and depth) and performing an affect on the DE variables upon interaction.

Note that it is easier to use the wrapper struct `Ray` to compute the solution, instead of calling this function.
"""
function propagation_problem(
	θ₀::Real,
	src::Source,
	ocn::Medium,
	bty::Boundary,
	ati::Boundary)

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

		∂²c_∂n²(r, z) = ocn.c(r, z)^2*(
			ocn.∂²c_∂r²(r, z)*ζ^2
			- 2ocn.∂²c_∂r∂z(r, z)*ξ*ζ
			+ ocn.∂²c_∂z²(r, z)*ξ^2
		)

		du[1] = dr_ds = ocn.c(r, z)*ξ
		du[2] = dz_ds = ocn.c(r, z)*ζ
		du[3] = dξ_ds = -ocn.∂c_∂r(r, z)/ocn.c(r, z)^2
		du[4] = dζ_ds = -ocn.∂c_∂z(r, z)/ocn.c(r, z)^2
		du[5] = dτ_ds = 1/ocn.c(r, z)
		du[6] = dpʳ_ds = ∂²c_∂n²(r, z)/ocn.c(r, z)^2*qʳ
		du[7] = dpⁱ_ds = ∂²c_∂n²(r, z)/ocn.c(r, z)^2*qⁱ
		du[8] = dqʳ_ds = ocn.c(r, z)*pʳ
		du[9] = dqⁱ_ds = ocn.c(r, z)*pⁱ
	end

	rng_condition(u, t, ray) = ocn.R/2 - abs(u[1] - ocn.R/2)
	rng_affect!(ray) = terminate!(ray)
	CbRng = ContinuousCallback(rng_condition, rng_affect!)
	CbBty = ContinuousCallback(bty.condition, bty.affect!)
	CbAti = ContinuousCallback(ati.condition, ati.affect!)
	CbBnd = CallbackSet(CbRng, CbBty, CbAti)

	r₀ = src.pos.r
	z₀ = src.pos.z
	ξ₀ = cos(θ₀)/ocn.c(r₀, z₀)
	ζ₀ = sin(θ₀)/ocn.c(r₀, z₀)
	τ₀ = 0.0

	λ₀ = ocn.c(r₀, z₀)/src.sig.f
	ω = src.sig.f
	p₀ʳ = 1.0
	p₀ⁱ = 0.0
	W₀ = 100λ₀ # 10..50
	q₀ʳ = 0.0
	q₀ⁱ = ω*W₀^2/2

	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀ʳ, p₀ⁱ, q₀ʳ, q₀ⁱ]

	TLmax = 100.0
	S = 10^(TLmax/10.0)
	sSpan = (0., S)

	prob = ODEProblem(propagation!, u₀, sSpan)

	return prob, CbBnd
end

"""
	solve_propagation(prob_eikonal::ODEProblem, CbBnd::ContinuousCallback) -> RaySol::ODESolution

Solves the eikonal, transport, and dynamic ray equations defined by the differential equation problem `prob` with continuous callback `CbBnd`.

Returns the `ODESolution` `RaySol` which is a type belonging to the `DifferentialEquations.jl` package.

Note that it is easier to use the wrapper struct `Ray` to compute the solution, instead of calling this function.
"""
function solve_propagation(prob::ODEProblem, CbBnd::Union{ContinuousCallback, CallbackSet})
	RaySol = solve(prob, AutoVern7(Rodas4()), callback = CbBnd, reltol=1e-8, abstol=1e-8)
	return RaySol
end

struct Ray
	θ₀::Real
	sol
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
end

"""
	Ray(θ₀::Real, src::Source, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0.0))

Computes a ray path launched from `src` at an initial angle `θ₀` within an ocean medium `ocn` bounded above by the altimetry `ati` and bathymetry `bty`. The default altimetry is a flat sea surface.

The following fields are stored in an instance of `Ray`:
* `θ₀` the initial ray angle (radians)
* `sol` the `ODESolution`
* `S` the maximum ray path `s` length (metres)
* `r(s)` range (metres) as a function of arc length `s` (metres)
* `z(s)` depth (metres) as a function of arc length `s` (metres)
* `ξ(s)` range component ray slowness (s/m) as a function of arc length `s` (metres)
* `ζ(s)` depth component ray slowness (s/m) as a function of arc length `s` (metres)
* `τ(s)` time lapsed of ray journey as a function of arc length `s` (metres)
* `p(s)` slowness (s/m) as a function of arc length `s` (metres)
* `q(s)` spreading (m/rad) as a function of arc length `s` (metres)
* `θ(s)` as a function of arc length `s` (metres)
* `c(s)` as a function of arc length `s` (metres)
"""
function Ray(θ₀::Real, src::Source, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0.0))
	Prob, CbBnd = propagation_problem(θ₀, src, ocn, bty, ati)
	sol = solve_propagation(Prob, CbBnd)

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

	return Ray(θ₀, sol, S, r, z, ξ, ζ, τ, p, q, θ, c)
end

struct Beam
	ray
	b::Function
	W::Function
end

"""
	Beam(θ₀::Real, src::Source, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0.0))

Computes a complex-valued Gaussian pressure beam propagating through space for the `Ray` trace solved by the scenario defined by the input parameters.

The fields stored are:
* `θ₀` initial ray angle (radians)
* `ray` the `Ray` solution struct
* `b(s, n)` the complex-valued pressure beam (Pa) as a function of arc length `s` (metres) and arc normal `n` (metres)
* `S` maximum arc length (metres)
* `W(s)` the computed half-beamwidth (metres) in terms of arc length `s` (metres)
"""
function Beam(θ₀::Real, src::Source, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0.0))
	
	ray = Ray(θ₀, src, ocn, bty, ati)
	
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

	A = 1/c₀ * exp(im*π/4)*sqrt(q₀*ω*cos(θ₀)/2π)
	b(s, n) = A * sqrt(c(s)/r(s)/q(s)) * exp(-im*ω * (τ(s) + p(s)/q(s)*n^2/2))

	return Beam(ray, b, W)
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

Base.broadcastable(m::Position) = Ref(m)
Base.broadcastable(m::Medium) = Ref(m)
Base.broadcastable(m::Boundary) = Ref(m)
Base.broadcastable(m::Signal) = Ref(m)
Base.broadcastable(m::Source) = Ref(m)
Base.broadcastable(m::Ray) = Ref(m)

function add_to_pressure(r::Real, z::Real, beam::Beam, δθ₀::Real, coh_pre::Function)
	sMins, nMins = closest_points(r, z, beam)
	p = complex(0)
	for (n, sMin) ∈ enumerate(sMins)
		p += coh_pre(δθ₀ * beam.b(sMin, nMins[n]))
	end
	return p
end

struct Field
	beams::AbstractVector{Beam}
	src::Source
	ocn::Medium
	bty::Boundary
	ati::Boundary
	p::Function
	TL::Function
end

function Field(beams::AbstractVector{T}, src::Source, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0.0)) where T <: Beam
	δθ₀ = diff([beams[nBeam].ray.θ₀ for nBeam = eachindex(beams)])
	push!(δθ₀, δθ₀[end])

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

## Plots
function acoustic_plot()
	f = Figure()
	hold!(f, true)
	xlabel!(f, "Range (m)")
	ylabel!(f, "Depth (m)")
	colorscheme!(f, "light")
	return f
end

function acoustic_plot!(title::AbstractString)
	f = gcf()
	title!(f, title)
end

function acoustic_plot!()
	f = gcf()
	yflip!(f, true)
end

function acoustic_plot!(bnd::Boundary)
	f = gcf()
	r = LinRange(0.0, bnd.R, 1001)
	plot!(f, r, bnd.z,
		linecolor = color(0, 0, 0))
end

function acoustic_plot(bnd::Boundary)
	f = acoustic_plot()
	acoustic_plot!(bnd)
	acoustic_plot!()
	return f
end

function acoustic_plot!(ray::Ray)
	f = gcf()
	s = LinRange(0.0, ray.S, 1001)
	plot!(f, ray.r(s), ray.z(s))
end

function acoustic_plot(ray::Ray)
	f = acoustic_plot()
	acoustic_plot!(ray)
	acoustic_plot!()
	return f
end

function acoustic_plot!(fld::Field)
	f = gcf()
	r = LinRange(0.0, fld.ocn.R, 11)
	z = LinRange(0.0, fld.ocn.Z, 5)
	contourf(r, z, fld.TL.(r', z),
		levels = 21, majorlevels = 2)
end

function acoustic_plot(fld::Field)
	f = acoustic_plot()
	acoustic_plot!(fld)
	acoustic_plot!()
	return f
end
