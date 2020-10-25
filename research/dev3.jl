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
AutoVern7,
ODECompositeSolution
using ForwardDiff:
derivative,
gradient
using Base: broadcastable
using Roots: find_zeros
using IntervalArithmetic:
(..),
Interval

abstract type OceanAcoustic <: Any end

function Base.show(io::IO, oac::OceanAcoustic)
	println(io, string(typeof(oac)), "(")
	for p ∈ propertynames(oac)
		println(io, " ", p, "::", typeof(getproperty(oac, p)))
	end
	print(io, ")")
end

Base.broadcastable(oac::OceanAcoustic) = Ref(oac)

parse_interval(Ω::Interval) = Ω
parse_interval(V::Real) = V..V

function interpolated_function(x, y)
	Itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
	return ItpFcn(x::Real) = Itp(x)
end
function interpolated_function(x, y, z)
	Itp = LinearInterpolation((x, y), z, extrapolation_bc = Flat())
	return ItpFcn(x::Real, y::Real) = Itp(x, y)
end

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
	callback::ContinuousCallback
	R::Real
	
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
		callback = ContinuousCallback(condition, affect!)
		return new(z, dz_dr, callback, R)
	end
end

function Boundary(
	r::Vector{T},
	z::Vector{T},
	R::T = r[end]) where T <: Real
	zFcn = interpolated_function(r, z)
	return Boundary(zFcn, R)
end

function Boundary(rz::AbstractArray{T}, R::T = rz[end, 1]) where T <: Real
	r = [rng for rng ∈ rz[:, 1]]
	z = [dpt for dpt ∈ rz[:, 2]]
	return Boundary(r, z, R)
end

function Boundary(z::Real, R::Real = 1.0)
	zFcn(r) = z
	return Boundary(zFcn, R)
end

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

struct Medium <: OceanAcoustic
	SSP::Celerity
	R::Real
	Z::Real

	function Medium(c::Function, R::Real, Z::Real)
		SSP = Celerity(c)
		return new(SSP, R, Z)
	end
end

function Medium(c::AbstractArray, R::Real = c[end, 1], Z::Real = c[1, end])
	r_ = [rc for rc ∈ c[1, 2:end]]
	z_ = [zc for zc ∈ c[2:end, 1]]
	c_ = c[2:end, 2:end]'
	
	cFcn = interpolated_function(r_, z_, c_)
	return Medium(cFcn, R, Z)
end

function Medium(z::AbstractVector, c::AbstractVector, R::Real, Z::Real = z[end])
	cMat = vcat([0 0 R], hcat(z, c, c))
	return Medium(cMat, R, Z)
end

function Medium(c::Real, R::Real, Z::Real)
	cFcn(r, z) = c
	return Medium(cFcn, R, Z)
end

struct Environment <: OceanAcoustic
	Ωr::Interval
	Ωz::Interval
	callback::ContinuousCallback
	ocn::Medium
	bty::Boundary
	ati::Boundary
	function Environment(Ωr::Interval, ocn::Medium, bty::Boundary, ati::Boundary = Boundary(0))
		bty_Ωz = bty.z(Ωr) |> parse_interval
		ati_Ωz = ati.z(Ωr) |> parse_interval
		Ωz = bty_Ωz ∪ ati_Ωz

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

function propagation(scn::Scenario)
	c(r, z) = scn.env.ocn.SSP.c(r, z)
	∂c_∂r(r, z) = scn.env.ocn.SSP.∂c_∂r(r, z)
	∂c_∂z(r, z) = scn.env.ocn.SSP.∂c_∂z(r, z)
	∂²c_∂r²(r, z) = scn.env.ocn.SSP.∂²c_∂r²(r, z)
	∂²c_∂z²(r, z) = scn.env.ocn.SSP.∂²c_∂z²(r, z)
	∂²c_∂r∂z(r, z) = scn.env.ocn.SSP.∂²c_∂r∂z(r, z)

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
	for θ₀ ∈ scn.src.fan.θ₀s
		ξ₀ = cos(θ₀) / c(r₀, z₀)
		ζ₀ = sin(θ₀) / c(r₀, z₀)
		u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, pʳ₀, pⁱ₀, qʳ₀, qⁱ₀]
	
		prob = ODEProblem(propagation!, u₀, sSpan)
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

	return sols
end

struct Ray <: OceanAcoustic
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

function convergence()
	c = [1520, 1500, 1515, 1495, 1545.]
	z = [0., 300., 1200., 2e3, 5000.]
	Z = z[end]
	R = 250e3

	ocn = Medium(z, c, R, Z)
	bty = Boundary(5e3, R)
	ati = Boundary(0., R)
	env = Environment(R, ocn, bty, ati)

	θ_crit = acos(ocn.SSP.c(0.0, 0.0)/ocn.SSP.c(0.0, 5e3))
	θ₀s = θ_crit*LinRange(0.5, 1.0, 10)

	fan = Fan(θ₀s)
	src = Source(Position(0., 0.), Signal(200.), fan)
	scn = Scenario(env, src, "Convergence Zone Propagation")
end

function wavy()
	# Altimetry
	zAtiMin = 0.
	zAtiMax = 50.
	zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2

	# Bathymetry
	rBtyPeak = 5e3
	zBtyMax = 1e3
	zBtyMin = 8e2
	Aᵣ = (2rBtyPeak/3)^2/log((9.0zBtyMax - 11.0zBtyMin)/(10.0(zBtyMax - zBtyMin)))
	zBty(r) = zBtyMax - (zBtyMax - zBtyMin)*exp(-(r - rBtyPeak)^2/4e5)

	# Ocean
	rOcnMax = 10e3
	cOcnMin = 1500.
	cOcnMax = 1600.
	cSolve(r) = [
		1.0 zAti(r) zAti(r)^2
		1.0 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2)^2
		1.0 zBty(r) zBty(r)^2
	]
	cCoeff(r) = cSolve(r)\[cOcnMax, cOcnMin, cOcnMax]
	# cCoeff₀(r) = cCoeff(r)[1]
	# cCoeff₁(r) = cCoeff(r)[2]
	# cCoeff₂(r) = cCoeff(r)[3]
	# cOcn(r, z) = cCoeff₀(r) + cCoeff₁(r)*z + cCoeff₂(r)*z^2
	cOcn(r, z) = cCoeff(r)' * [1, z, z^2]

	# Environment
	ocn = Medium(cOcn, rOcnMax, zBtyMax)
	bty = Boundary(zBty, rOcnMax)
	ati = Boundary(zAti, rOcnMax)
	env = Environment(rOcnMax, ocn, bty, ati)

	# Scenario
	r₀ = 0.0
	z₀ = (zBty(r₀) + zAti(r₀))/2
	θ_crit = acos(cOcn(r₀, z₀)/cOcnMax)

	fan = Fan(θ_crit*(-1.5:0.125:1.5))
	src = Source(Position(r₀, z₀), Signal(250.), fan)
	scn = Scenario(env, src, "Wavy Environment")
end

##
scn = wavy()

sols = propagation(scn)

##
using Plots

p = plot(yaxis = :flip)
plot!.(sols, vars = (1, 2))
display(p)