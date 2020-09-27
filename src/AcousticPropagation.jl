import Interpolations: LinearInterpolation
import DifferentialEquations:
ContinuousCallback,
CallbackSet,
ODEProblem,
solve,
terminate!
import ForwardDiff: ForwardDiff
import Base: broadcastable

export Position
export Signal
export Source
export Boundary
export Medium
export Ray
export Beam

function InterpolatingFunction(rng, val)
	Itp = LinearInterpolation(rng, val, extrapolation_bc = Flat())
	return ItpFcn(r) = Itp(r)
end
function InterpolatingFunction(rng, dpt, val)
	Itp = LinearInterpolation((dpt, rng), val, extrapolation_bc = Flat())
	return ItpFcn(r, z) = Itp(z, r)
end

function boundary_reflection(t_inc::Vector, t_bnd::Vector)
	# works for parabolic boundary
	MyAngle(tng) = atan(tng[2]/tng[1])
	θ_inc = MyAngle(t_inc)
	θ_bnd = MyAngle(t_bnd)

	c = cos(θ_inc)/t_inc[1]

	θ_inc_flat = θ_inc - θ_bnd
	θ_rfl_flat = -θ_inc_flat
	θ_rfl = θ_rfl_flat + θ_bnd

	return [cos(θ_rfl), sin(θ_rfl)]/c
end

struct Position
	r::Real
	z::Real
end

struct Signal
	f::Real
end

struct Source
	Pos::Position
	Sig::Signal
end

struct Boundary
	z::Function
	dz_dr::Function
	condition::Function
	affect!::Function
	function Boundary(z::Function)
		dz_dr(r) = ForwardDiff.derivative(z, r)
		condition(u, t, ray) = z(u[1]) - u[2]
		affect!(ray) = ray.u[3], ray.u[4] = boundary_reflection([ray.u[3], ray.u[4]], [1, dz_dr(ray.u[1])])
		return new(z, dz_dr, condition, affect!)
	end
end
function Boundary(r::Vector, z::Vector)
	zFcn = InterpolatingFunction(r, z)
	return Boundary(zFcn)
end
function Boundary(rz::AbstractArray)
	r = [rng for rng ∈ rz[:, 1]]
	z = [dpt for dpt ∈ rz[:, 2]]
	return Boundary(r, z)
end
function Boundary(z::Real)
	zFcn(r) = z
	return Boundary(zFcn)
end

struct Medium
	c::Function
	∂c_∂r::Function
	∂c_∂z::Function
	∂²c_∂r²::Function
	∂²c_∂r∂z::Function
	∂²c_∂z²::Function
	R::Real
	function Medium(c::Function, R::Real)
		c_(x) = c(x[1], x[2])
		∇c_(x) = ForwardDiff.gradient(c_, x)
		∇c(r, z) = ∇c_([r, z])
		∂c_∂r(r, z) = ∇c(r, z)[1]
		∂c_∂z(r, z) = ∇c(r, z)[2]
	
		∂c_∂r_(x) = ∂c_∂r(x[1], x[2])
		∇∂c_∂r_(x) = ForwardDiff.gradient(∂c_∂r_, x)
		∇∂c_∂r(r, z) = ∇∂c_∂r_([r, z])
	
		∂c_∂z_(x) = ∂c_∂z(x[1], x[2])
		∇∂c_∂z_(x) = ForwardDiff.gradient(∂c_∂r_, x)
		∇∂c_∂z(r, z) = ∇∂c_∂z_([r, z])
	
		∂²c_∂r²(r, z) = ∇∂c_∂r(r, z)[1]
		∂²c_∂r∂z(r, z) = ∇∂c_∂r(r, z)[2]
		∂²c_∂z²(r, z) = ∇∂c_∂z(r, z)[2]
	
		return new(c, ∂c_∂r, ∂c_∂z, ∂²c_∂r², ∂²c_∂r∂z, ∂²c_∂z², R)
	end
end
function Medium(c::AbstractArray, R::Real = c[end, 1])
	r_ = [rc for rc ∈ c[1, 2:end]]
	z_ = [zc for zc ∈ c[2:end, 1]]
	c_ = c[2:end, 2:end]
	
	cFcn = InterpolatingFunction(r_, z_, c_)
	return Medium(cFcn, R)
end
function Medium(z::AbstractVector, c::AbstractVector, R = z[end])
	cMat = vcat([0 0 R], hcat(z, c, c))
	return Medium(cMat, R)
end
function Medium(c::Real, R::Real)
	cFcn(r, z) = c
	return Medium(cFcn, R)
end

function acoustic_propagation_problem(
	θ₀::Real,
	Src::Source,
	Ocn::Medium,
	Bty::Boundary,
	Ati::Boundary)

	function eikonal!(du, u, p, s)
		r = u[1]
		z = u[2]
		ξ = u[3]
		ζ = u[4]
		τ = u[5]
		pʳ = u[6]
		pⁱ = u[7]
		qʳ = u[8]
		qⁱ = u[9]

		∂²c_∂n²(r, z) = Ocn.c(r, z)^2*(
			Ocn.∂²c_∂r²(r, z)*ζ^2
			- 2Ocn.∂²c_∂r∂z(r, z)*ξ*ζ
			+ Ocn.∂²c_∂z²(r, z)*ξ^2
		)

		du[1] = dr_ds = Ocn.c(r, z)*ξ
		du[2] = dz_ds = Ocn.c(r, z)*ζ
		du[3] = dξ_ds = -Ocn.∂c_∂r(r, z)/Ocn.c(r, z)^2
		du[4] = dζ_ds = -Ocn.∂c_∂z(r, z)/Ocn.c(r, z)^2
		du[5] = dτ_ds = 1/Ocn.c(r, z)
		du[6] = dpʳ_ds = ∂²c_∂n²(r, z)/Ocn.c(r, z)^2*qʳ
		du[7] = dpⁱ_ds = ∂²c_∂n²(r, z)/Ocn.c(r, z)^2*qⁱ
		du[8] = dqʳ_ds = Ocn.c(r, z)*pʳ
		du[9] = dqⁱ_ds = Ocn.c(r, z)*pⁱ
	end

	rng_condition(u, t, ray) = Ocn.R/2 - abs(u[1] - Ocn.R/2)
	rng_affect!(ray) = terminate!(ray)
	CbRng = ContinuousCallback(rng_condition, rng_affect!)
	CbBty = ContinuousCallback(Bty.condition, Bty.affect!)
	CbAti = ContinuousCallback(Ati.condition, Ati.affect!)
	CbBnd = CallbackSet(CbRng, CbBty, CbAti)

	r₀ = Src.Pos.r
	z₀ = Src.Pos.z
	ξ₀ = cos(θ₀)/Ocn.c(r₀, z₀)
	ζ₀ = sin(θ₀)/Ocn.c(r₀, z₀)
	τ₀ = 0.0

	λ₀ = Ocn.c(r₀, z₀)/Src.Sig.f
	ω = Src.Sig.f
	p₀ʳ = 1.0
	p₀ⁱ = 0.0
	W₀ = 100λ₀ # 10..50
	q₀ʳ = 0.0
	q₀ⁱ = ω*W₀^2/2

	u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀ʳ, p₀ⁱ, q₀ʳ, q₀ⁱ]

	TLmax = 100
	S = 10^(TLmax/10)
	sSpan = (0., S)

	prob_eikonal = ODEProblem(eikonal!, u₀, sSpan)

	return prob_eikonal, CbBnd
end

function solve_acoustic_propagation(prob_eikonal, CbBnd)
	@time RaySol = solve(prob_eikonal, callback = CbBnd, reltol=1e-8, abstol=1e-8)
	return RaySol
end

struct Ray
	θ₀
	Sol
	S
	r
	z
	ξ
	ζ
	τ
	p
	q
	θ
	c
	function Ray(θ₀::Real, Src::Source, Ocn::Medium, Bty::Boundary, Ati::Boundary = Boundary(0))
		Prob, CbBnd = acoustic_propagation_problem(θ₀, Src, Ocn, Bty, Ati)
		Sol = solve_acoustic_propagation(Prob, CbBnd)
	
		S = Sol.t[end]
		r(s) = Sol(s, idxs = 1)
		z(s) = Sol(s, idxs = 2)
		ξ(s) = Sol(s, idxs = 3)
		ζ(s) = Sol(s, idxs = 4)
		τ(s) = Sol(s, idxs = 5)
		p(s) = Sol(s, idxs = 6) + im*Sol(s, idxs = 7)
		q(s) = Sol(s, idxs = 8) + im*Sol(s, idxs = 9)
		θ(s) = atan(ζ(s)/ξ(s))
		c(s) = cos(θ(s))/ξ(s)
	
		return new(θ₀, Sol, S, r, z, ξ, ζ, τ, p, q, θ, c)
	end
end

struct Beam
	θ₀::Real
	ray
	b::Function
	S::Real
	W::Function
	function Beam(θ₀::Real, Src::Source, Ocn::Medium, Bty::Boundary, Ati::Boundary = Boundary(0))
		
		ray = Ray(θ₀, Src, Ocn, Bty, Ati)
		
		r(s) = ray.r(s)
		z(s) = ray.z(s)
		τ(s) = ray.τ(s)
		p(s) = ray.p(s)
		q(s) = ray.q(s)
		c(s) = ray.c(s)
		W(s) = sqrt(-2/ω/imag(p(s)/q(s)))
	
		c₀ = c(0)
		ω = 2π*Src.Sig.f
		λ₀ = c₀/Src.Sig.f
		W₀ = W(0)
		q₀ = q(0)
	
		A = 1/c₀ * exp(im*π/4)*sqrt(q₀*ω*cos(θ₀)/2π)
		b(s, n) = A * sqrt(c(s)/r(s)/q(s)) * exp(-im*ω * (τ(s) + p(s)/q(s)*n^2/2))

		return new(θ₀, ray, b, ray.S, W)
	end
end

# function addtofield!(p, r, z, b)

# end

# struct Field
# 	θ₀::Union{Real,Vector}
# 	p::Func
# 	function Field(
# 		θ₀vals::Vector,
# 		Src::Source,
# 		Rcv::Receiver,
# 		Ocn::Medium,
# 		Bty::Boundary,
# 		Ati::Boundary = Boundary(0),
# 		Before::Function = p -> p,
# 		After!::Function = p -> p)

# 		p = zeros(length(Rcv.r), length(Rcv.z))
# 		θ₀s = sort(θ₀vals)
# 		rays = Beam.(θ₀s, Src, Ocn, Bty, Ati)
		
# 		θ₀ = θ₀s[1]
# 		δθ = θ₀ - θ₀s[2]
# 		b(s, n) = Before(δθ*rays[1].b(s, n))
# 		addtofield!(p, Rcv.r, Rcv.z, b)

# 		for n = 2:length(θ₀)-1

# 		end

# 		After!(p)
# 		return new(θ₀, p, TL)
# 	end
# end

Base.broadcastable(m::Position) = Ref(m)
Base.broadcastable(m::Medium) = Ref(m)
Base.broadcastable(m::Boundary) = Ref(m)
Base.broadcastable(m::Signal) = Ref(m)
Base.broadcastable(m::Source) = Ref(m)
# Base.broadcastable(m::Ray) = Ref(m)
