## Preamble
using Interpolations:
LinearInterpolation,
Flat
using IntervalArithmetic
using ForwardDiff: derivative
using DifferentialEquations
using Plots

## Auxiliaries
parse_interval(Ω::Interval) = Ω
parse_interval(V::Real) = V..V

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

## Environment
# Slice
R = 250e3

# Medium: Ocean
cVec = [1520, 1500, 1515, 1495, 1545]
zcVec = [0, 300, 1200, 2e3, 5000]

cItp = LinearInterpolation(zcVec, cVec, extrapolation_bc = Flat())
c(r, z) = cItp(z)

∂c_∂r(r, z) = derivative(r -> c(r, z), r)
∂c_∂z(r, z) = derivative(z -> c(r, z), z)

∂²c_∂r²(r, z) = derivative(r -> ∂c_∂r(r, z), r)
∂²c_∂z²(r, z) = derivative(z -> ∂c_∂z(r, z), z)
∂²c_∂r∂z(r, z) = derivative(r -> ∂c_∂z(r, z), r)
∂²c_∂z∂r(r, z) = derivative(z -> ∂c_∂r(r, z), z)

# Boundary: Bathymetry
zBty(r) = 5000.
dzBty_dr(r) = derivative(zBty, r)
bty_condition(u, t, ray) = zBty(u[1]) - u[2]
function bty_affect!(ray)
	rᵢ = ray.u[1]
	ξᵢ = ray.u[3]
	ζᵢ = ray.u[4]

	ξₒ, ζₒ = boundary_reflection(
		[ξᵢ, ζᵢ],
		[1, dzBty_dr(rᵢ)]
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
bty_callback = ContinuousCallback(bty_condition, bty_affect!)

# Boundary: Altimetry
zAti(r) = 0.
dzAti_dr(r) = derivative(zAti, r)
ati_condition(u, t, ray) = zAti(u[1]) - u[2]
function ati_affect!(ray)
	rᵢ = ray.u[1]
	ξᵢ = ray.u[3]
	ζᵢ = ray.u[4]

	ξₒ, ζₒ = boundary_reflection(
		[ξᵢ, ζᵢ],
		[1, dzAti_dr(rᵢ)]
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
ati_callback = ContinuousCallback(ati_condition, ati_affect!)

# Slice
Ωr = 0..R
Ωz_bty = Ωr |> zBty |> parse_interval
Ωz_ati = Ωr |> zAti |> parse_interval
Ωz = Ωz_bty ∪ Ωz_ati

env_condition(u, t, ray) = (Ωr.hi - Ωr.lo)/2 - abs(u[1] - (Ωr.hi + Ωr.lo)/2)
env_affect!(ray) = terminate!(ray)
env_callback = ContinuousCallback(env_condition, env_affect!)

## Scenario
# Source
r₀ = 0.
z₀ = 0.
f = 200.
ω = 2π*f

# Spark
θ_crit = c(r₀, z₀)/c(r₀, Ωz.hi) |> acos

# Fan
θ₀s = θ_crit * LinRange(0, 1, 11)
δθ₀s = diff(θ₀s)

## Trace
callbacks = CallbackSet(env_callback, bty_callback, ati_callback)

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

	∂²c_∂n²(r, z) = c(r, z)^2*(
		∂²c_∂r²(r, z)*ζ^2
		- 2∂²c_∂r∂z(r, z)*ξ*ζ
		+ ∂²c_∂z²(r, z)*ξ^2
	)

	du[1] = dr_ds = c(r, z)*ξ
	du[2] = dz_ds = c(r, z)*ζ
	du[3] = dξ_ds = -∂c_∂r(r, z)/c(r, z)^2
	du[4] = dζ_ds = -∂c_∂z(r, z)/c(r, z)^2
	du[5] = dτ_ds = 1/c(r, z)
	du[6] = dpʳ_ds = ∂²c_∂n²(r, z)/c(r, z)^2*qʳ
	du[7] = dpⁱ_ds = ∂²c_∂n²(r, z)/c(r, z)^2*qⁱ
	du[8] = dqʳ_ds = c(r, z)*pʳ
	du[9] = dqⁱ_ds = c(r, z)*pⁱ
end

sols = Vector{DiffEqBase.AbstractODESolution}(undef, 0)
for θ₀ ∈ θ₀s
	ξ₀ = cos(θ₀)/c(r₀, z₀)
	ζ₀ = sin(θ₀)/c(r₀, z₀)
	τ₀ = 0.0

	λ₀ = c(r₀, z₀)/f
	p₀ʳ = 1.0
	p₀ⁱ = 0.0
	W₀ = 100λ₀ # 10..50
	q₀ʳ = 0.0
	q₀ⁱ = ω*W₀^2/2

	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀ʳ, p₀ⁱ, q₀ʳ, q₀ⁱ]

	TLmax = 100.
	S = 10^(TLmax/10)
	sSpan = (0., S)

	prob = ODEProblem(propagation!, u₀, sSpan)
	push!(sols, solve(prob, callback = callbacks))
end

## Plot
using Plots

r = LinRange(Ωr.lo, Ωr.hi, 101)
z = LinRange(Ωz.lo, Ωz.hi, 51)

p = plot()
plot!(r, zBty)
plot!(r, zAti)
plot!.(sols, vars = (1, 2))
display(p)
