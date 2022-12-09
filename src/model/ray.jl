dot(u::AbstractVector{<:Number}, v::AbstractVector{<:Number}) = u' * v

function gradient2tangent(m::Float64)
	x = 1 / √(m^2 + 1)
	[x; m*x]
end

function reflection(tng_inc::Vector{Float64}, m_bnd::Float64)
	tng_bnd = gradient2tangent(m_bnd)
	nrm_bnd = [0 1; -1 0] * tng_bnd
	α = dot(tng_inc, nrm_bnd)
	tng_rfl = tng_inc - 2α * nrm_bnd
	atan(tng_rfl[2] / tng_rfl[1])
end

"""
`Ray`

Fields:
* `θ₀::Float64` initial ray angle
* `s_max::Float64` maximum ray length
* `r::Function` Univariate range of ray wrt ray length
* `z::Function` Univariate depth of ray wrt ray length
* `p::Function`:
	* Univariate ray pressure `p(r)`
	* Bivariate beam pressure `p(r, z)`
"""
struct Ray <: Oac
	θ₀::Float64
	s_max::Float64
	r::Function
	z::Function
	c::Function
	τ::Function
	p::Function
	TL::Function
end

export Ray

# struct Rays <: Oac
# 	vec_of_rays::Vector{Ray}
# 	function Rays(rays::AbstractVector{Ray})
# 		new(rays)
# 	end
# end

# export Rays

# function Base.getindex(rays::Rays, n::Int)
# 	rays.vec_of_rays[n]
# end

# Base.firstindex(rays::Rays) = 1
# Base.lastindex(rays::Rays) = length(rays.vec_of_rays)

function default_angles(scn::Scenario, N::Int = 1001)
	x_rcv = scn.ent.rcv.x
	rng_ocn = calc_ocean_depth_range(scn)

	dz_ocn = diff(rng_ocn)[1]
	θ₀ = atan(dz_ocn / (x_rcv / 10))

	return if scn.ent.src.z == scn.env.srf.z(0.0)
		θ₀ * range(0, 1, N)
	elseif scn.ent.src.z == scn.env.btm.z(0.0)
		θ₀ * range(-1, 0, N)
	else
		θ₀ * range(-1, 1, N)
	end
end

# struct Trace <: Oac
# 	scn::Scenario
# 	rays::Vector{Ray}
# end

struct Trace <: Oac
	rays::Vector{Ray}
end

export Trace

# function Trace(
# 	scn::Scenario,
# 	angles::AbstractVector{<:AbstractFloat} = default_angles(scn)
# 	)
# 	# Computation

# 	cel(r, z) = scn.env.ocn.c(r, z)
# 	∂c_∂r(r, z) = derivative(r̂ -> cel(r̂, z), r)
# 	∂c_∂z(r, z) = derivative(ẑ -> cel(r, ẑ), z)
# 	∂²c_∂r²(r, z) = derivative(r̂ -> ∂c_∂r(r̂, z), r)
# 	∂²c_∂z²(r, z) = derivative(ẑ -> ∂c_∂z(r, ẑ), z)
# 	∂²c_∂r∂z(r, z) = derivative(r̂ -> ∂c_∂z(r̂, z), r)
# 	∂²c_∂n²(r, z, ξ, ζ) = cel(r, z)^2 * (
# 		∂²c_∂r²(r, z) * ζ^2
# 		 - 2 * ∂²c_∂r∂z(r, z) * ξ * ζ
# 		 + ∂²c_∂z²(r, z) * ξ^2
# 	)

# 	δθ₀ = if length(angles) == 1
# 		1.0
# 	else
# 		angles |> diff |> mean
# 	end

# 	function reflect!(int, bnd)
# 		r, z, ξ, ζ = int.u[1:4]
# 		c = cel(r, z)
# 		tng_inc = c * [ξ; ζ]

# 		θ_rfl = reflection(tng_inc, derivative(bnd.z, r))

# 		int.u[3] = cos(θ_rfl) / c
# 		int.u[4] = sin(θ_rfl) / c
# 	end

# 	cb_rng = ContinuousCallback(
# 		(x, s, int) -> x[1] - scn.ent.rcv.x,
# 		terminate!
# 	)

# 	cb_srf = ContinuousCallback(
# 		(x, s, int) -> x[2] - scn.env.srf.z(x[1]),
# 		int -> reflect!(int, scn.env.srf)
# 	)

# 	cb_btm = ContinuousCallback(
# 		(x, s, int) -> x[2] - scn.env.btm.z(x[1]),
# 		int -> reflect!(int, scn.env.btm)
# 	)

# 	cb = CallbackSet(cb_rng, cb_btm, cb_srf)

# 	function tracer!(du, u, params, s)
# 		r, z, ξ, ζ, τ, pRe, pIm, qRe, qIm = u
		
# 		c = cel(r, z)
# 		∂cr = ∂c_∂r(r, z)
# 		∂cz = ∂c_∂z(r, z)
# 		∂²cnn = ∂²c_∂n²(r, z, ξ, ζ)

# 		du[1] = dr = c * ξ
# 		du[2] = dz = c * ζ
# 		du[3] = dξ = -∂cr / c^2
# 		du[4] = dζ = -∂cz / c^2
# 		du[5] = dτ = 1/c
# 		du[6] = dpRe = -∂²cnn / c^2 * qRe
# 		du[7] = dpIm = -∂²cnn / c^2 * qIm
# 		du[8] = dqRe = c * pRe
# 		du[9] = dqIm = c * pIm
# 	end

# 	f = scn.ent.src.f

# 	r₀ = 0.0
# 	z₀ = scn.ent.src.z
# 	c₀ = cel(r₀, z₀)
# 	τ₀ = 0.0
# 	p₀Re = 1.0
# 	p₀Im = 0.0
# 	λ₀ = c₀ / f
# 	W₀ = 30λ₀
# 	q₀Re = 0.0
# 	q₀Im = pi * f * W₀^2

# 	rays = Ray[]
# 	for θ₀ in angles
# 		ξ₀ = cos(θ₀) / c₀
# 		ζ₀ = sin(θ₀) / c₀

# 		u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀Re, p₀Im, q₀Re, q₀Im]

# 		s_span = (0.0, 100e3)
	
# 		prob = ODEProblem(tracer!, u₀, s_span)
	
# 		sol = solve(prob, Tsit5(), callback = cb)

# 		s_max = sol.t[end]
# 		r(s) = sol(s, idxs = 1)
# 		z(s) = sol(s, idxs = 2)
# 		ξ(s) = sol(s, idxs = 3)
# 		ζ(s) = sol(s, idxs = 4)
# 		τ(s) = sol(s, idxs = 5)
# 		p(s) = sol(s, idxs = 6) + im * sol(s, idxs = 7)
# 		q(s) = sol(s, idxs = 8) + im * sol(s, idxs = 9)
		
# 		c(s) = cel(r(s), z(s))
# 		ω = 2π * f
# 		A = δθ₀ / c(0.0) * sqrt(
# 			q(0.0) * f * cos(θ₀)
# 		) * exp(im * π / 4)
# 		function p_beam(s::Real, n::Real)
# 			(s > s_max || s < 0) && return ComplexF64(0.0)
# 			r(s) < 0 && return ComplexF64(0.0)
# 			return A * sqrt(
# 					c(s) / r(s) / q(s)
# 				) * exp(
# 				-im * ω * (
# 					τ(s) + p(s) / 2 / q(s) * n^2
# 				)
# 			)
# 		end
# 		p_beam(s::Real) = p_beam(s, 0.0)

# 		if save_trace
# 			push!(rays, Ray(θ₀, s_max, r, z, c, τ, p_beam, TL))
# 		end
# 	end
# 	new(scn, rays)
# end

# function Trace(
# 	scn::Scenario,
# 	N::Int = 1001
# 	)
# 	Trace(scn, default_angles(scn, N))
# end

# User Recipe
@userplot RayTracePlot
@recipe function plot(rtp::RayTracePlot)
	# Parse Inputs
	trc = rtp.args[1]

	# Plot Design.
	legend --> :none

	yflip := true

	# Ray Trace
	for ray in trc.rays
		s = range(0.0, ray.s_max, 501)
		r = ray.r.(s)
		z = ray.z.(s)
		@series r, z
	end
end

export raytraceplot

# function FieldOld(
# 	trc::Trace,
# 	rGrid::AbstractVector{<:Real} = range(0.0, trc.scn.ent.rcv.x, 250),
# 	zGrid::AbstractVector{<:Real} = range(calc_ocean_depth_range(trc.scn)..., 150)
# 	)

# 	Nr = length(rGrid)
# 	Nz = length(zGrid)
# 	p = zeros(ComplexF64, Nr, Nz)
# 	Narc = max(101, floor(Nr / 3) |> Int)
# 	for ray in trc.rays
# 		arc = range(0.0, ray.s_max, Narc)
# 		for i = eachindex(arc)[begin+1:end]
# 			sᵢ₋₁, sᵢ = arc[i-1 : i]
# 			rᵢ₋₁, rᵢ = ray.r.([sᵢ₋₁, sᵢ])
# 			# zᵢ₋₁ = ray.z(sᵢ₋₁)
# 			zᵢ = ray.z(sᵢ)
# 			for (nr, r) in enumerate(rGrid)
# 				if !(rᵢ₋₁ .≤ r .< rᵢ)
# 					continue
# 				end
# 				for (nz, z) in enumerate(zGrid)
# 					x_rcv = [r, z]
# 					# x_ray = [rᵢ₋₁, zᵢ₋₁]
# 					x_ray = [rᵢ, zᵢ]
# 					t_ray = ray.c(sᵢ₋₁) * [ray.ξ(sᵢ₋₁), ray.ζ(sᵢ₋₁)]
# 					t_ray /= dot(t_ray, t_ray) |> sqrt
# 					n_ray = ray.c(sᵢ₋₁) * [-ray.ζ(sᵢ₋₁), ray.ξ(sᵢ₋₁)]
# 					n_ray /= dot(n_ray, n_ray) |> sqrt
# 					s = dot(x_rcv - x_ray, t_ray)
# 					n = dot(x_rcv - x_ray, n_ray) |> abs
# 					p_add = ray.p(sᵢ₋₁ + s, n)
# 					if !(isnan(p_add) || isinf(p_add))
# 						p[nr, nz] += p_add
# 					end
# 				end
# 			end
# 		end
# 	end
# 	TL = -20log10.(p .|> abs)
# 	TL = max.(TL, 0.0)
# 	TL = min.(TL, 100.0)
# 	rGrid, zGrid, TL
# end

# export FieldOld

default_ranges(scn::Scenario) = range(0.0, scn.ent.rcv.x, 250)

default_depths(scn::Scenario) = range(calc_ocean_depth_range(scn)..., 150)

function RayMethodField(scn::Scenario,
	angles::AbstractVector{<:AbstractFloat} = default_angles(scn);
	ranges::AbstractVector{<:AbstractFloat} = default_ranges(scn),
	depths::AbstractVector{<:AbstractFloat} = default_depths(scn),
	save_field::Bool = true,
	save_trace::Bool = false
	)

	if !(save_field || save_trace)
		error("No computation.") # flesh out
	end

	cel(r, z) = scn.env.ocn.c(r, z)
	∂c_∂r(r, z) = derivative(r̂ -> cel(r̂, z), r)
	∂c_∂z(r, z) = derivative(ẑ -> cel(r, ẑ), z)
	∂²c_∂r²(r, z) = derivative(r̂ -> ∂c_∂r(r̂, z), r)
	∂²c_∂z²(r, z) = derivative(ẑ -> ∂c_∂z(r, ẑ), z)
	∂²c_∂r∂z(r, z) = derivative(r̂ -> ∂c_∂z(r̂, z), r)
	∂²c_∂n²(r, z, ξ, ζ) = cel(r, z)^2 * (
		∂²c_∂r²(r, z) * ζ^2
		 - 2 * ∂²c_∂r∂z(r, z) * ξ * ζ
		 + ∂²c_∂z²(r, z) * ξ^2
	)

	δθ₀ = if length(angles) == 1
		1.0
	else
		angles |> diff |> mean
	end

	function reflect!(int, bnd)
		r, z, ξ, ζ = int.u[1:4]
		c = cel(r, z)
		tng_inc = c * [ξ; ζ]

		θ_rfl = reflection(tng_inc, derivative(bnd.z, r))

		int.u[3] = cos(θ_rfl) / c
		int.u[4] = sin(θ_rfl) / c
	end

	cb_rng = ContinuousCallback(
		(x, s, int) -> x[1] - scn.ent.rcv.x,
		terminate!
	)

	cb_srf = ContinuousCallback(
		(x, s, int) -> x[2] - scn.env.srf.z(x[1]),
		int -> reflect!(int, scn.env.srf)
	)

	cb_btm = ContinuousCallback(
		(x, s, int) -> x[2] - scn.env.btm.z(x[1]),
		int -> reflect!(int, scn.env.btm)
	)

	cb = CallbackSet(cb_rng, cb_btm, cb_srf)

	function tracer!(du, u, params, s)
		r, z, ξ, ζ, τ, pRe, pIm, qRe, qIm = u
		
		c = cel(r, z)
		∂cr = ∂c_∂r(r, z)
		∂cz = ∂c_∂z(r, z)
		∂²cnn = ∂²c_∂n²(r, z, ξ, ζ)

		du[1] = dr = c * ξ
		du[2] = dz = c * ζ
		du[3] = dξ = -∂cr / c^2
		du[4] = dζ = -∂cz / c^2
		du[5] = dτ = 1/c
		du[6] = dpRe = -∂²cnn / c^2 * qRe
		du[7] = dpIm = -∂²cnn / c^2 * qIm
		du[8] = dqRe = c * pRe
		du[9] = dqIm = c * pIm
	end

	f = scn.ent.src.f

	r₀ = 0.0
	z₀ = scn.ent.src.z
	c₀ = cel(r₀, z₀)
	τ₀ = 0.0
	p₀Re = 1.0
	p₀Im = 0.0
	λ₀ = c₀ / f
	W₀ = 30λ₀
	q₀Re = 0.0
	q₀Im = pi * f * W₀^2

	Nr = length(ranges)
	Nz = length(depths)
	pressure = zeros(ComplexF64, Nr, Nz)

	rays = Ray[]

	for θ₀ in angles
		ξ₀ = cos(θ₀) / c₀
		ζ₀ = sin(θ₀) / c₀

		u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀Re, p₀Im, q₀Re, q₀Im]

		s_span = (0.0, 100e3)
	
		prob = ODEProblem(tracer!, u₀, s_span)
	
		sol = solve(prob, Tsit5(), callback = cb)

		s_max = sol.t[end]
		r(s) = sol(s, idxs = 1)
		z(s) = sol(s, idxs = 2)
		ξ(s) = sol(s, idxs = 3)
		ζ(s) = sol(s, idxs = 4)
		τ(s) = sol(s, idxs = 5)
		p(s) = sol(s, idxs = 6) + im * sol(s, idxs = 7)
		q(s) = sol(s, idxs = 8) + im * sol(s, idxs = 9)
		
		c(s) = cel(r(s), z(s))
		ω = 2π * f
		A = δθ₀ / c(0.0) * sqrt(
			q(0.0) * f * cos(θ₀)
		) * exp(im * π / 4)
		function beam(s::Real, n::Real)
			(s > s_max || s < 0) && return ComplexF64(0.0)
			r(s) < 0 && return ComplexF64(0.0)
			return A * sqrt(
					c(s) / r(s) / q(s)
				) * exp(
				-im * ω * (
					τ(s) + p(s) / 2 / q(s) * n^2
				)
			)
		end
		beam(s::Real) = beam(s, 0.0)
		TL(s::Real, n::Real) = -20log10(beam(s, n) |> abs)
		TL(s::Real) = TL(s, 0.0)

		if save_trace
			push!(rays, Ray(θ₀, s_max, r, z, c, τ, beam, TL))
		end

		if save_field
			arc = range(0.0, s_max, max(101, floor(Nr/3) |> Int))
			for i = eachindex(arc)[begin+1 : end]
				sᵢ₋₁, sᵢ = arc[i-1 : i]
				rᵢ₋₁, rᵢ = r.([sᵢ₋₁, sᵢ])
				# zᵢ₋₁ = z(sᵢ₋₁)
				zᵢ = z(sᵢ)
				for (nr, r_grid) in enumerate(ranges)
					if !(rᵢ₋₁ .≤ r_grid .< rᵢ)
						continue
					end
					for (nz, z_grid) in enumerate(depths)
						x_rcv = [r_grid, z_grid]
						# x_ray = [rᵢ₋₁, zᵢ₋₁]
						x_ray = [rᵢ, zᵢ]
						t_ray = c(sᵢ₋₁) * [ξ(sᵢ₋₁), ζ(sᵢ₋₁)]
						t_ray /= dot(t_ray, t_ray) |> sqrt
						n_ray = c(sᵢ₋₁) * [-ζ(sᵢ₋₁), ξ(sᵢ₋₁)]
						n_ray /= dot(n_ray, n_ray) |> sqrt
						s = dot(x_rcv - x_ray, t_ray)
						n = dot(x_rcv - x_ray, n_ray) |> abs
						p_add = beam(sᵢ₋₁ + s, n)
						if !(isnan(p_add) || isinf(p_add))
							pressure[nr, nz] += p_add
						end
					end
				end
			end
		end
	end
	TL = -20log10.(pressure .|> abs)
	TL = max.(0.0, TL)
	TL = min.(100.0, TL)

	fld = Field(ranges, depths, TL)
	trc = Trace(rays)

	return if save_trace && save_field
		fld, trc
	elseif save_field
		fld
	elseif save_trace
		trc
	end
end

export RayMethod

function RayMethodField(scn::Scenario,
	Nangles::Int = 101;
	ranges::AbstractVector{<:Float64} = default_ranges(scn),
	depths::AbstractVector{<:Float64} = default_depths(scn),
	save_field::Bool = true,
	save_trace::Bool = false
	)

	return RayMethodField(scn,
		default_angles(scn, Nangles),
		ranges = ranges,
		depths = depths,
		save_field = save_field,
		save_trace = save_trace
	)
end

export RayMethodField
