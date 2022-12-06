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

struct Ray <: OAC
	θ₀
	r
	z
	s_max
end

function default_angles(scn::Scenario)
	x_rcv = scn.ent.rcv.x
	rng_ocn = [
		calc_bnd_range(scn, :srf).lo
		calc_bnd_range(scn, :btm).hi
	]

	dz_ocn = diff(rng_ocn)[1]
	θ₀ = atan(dz_ocn / (x_rcv / 10))

	return if scn.ent.src.z == scn.env.srf.z(0.0)
		θ₀ * range(0, 1, 21)
	elseif scn.ent.src.z == scn.env.btm.z(0.0)
		θ₀ * range(-1, 0, 21)
	else
		θ₀ * range(-1, 1, 21)
	end
end

struct Trace <: OAC
	scn::Scenario
	rays::Vector{Ray}

	function Trace(
		scn::Scenario,
		angles::AbstractVector{<:AbstractFloat} = default_angles(scn)
		)
		# Computation
	
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
		W₀ = 25λ₀
		q₀Re = 0.0
		q₀Im = pi * f * W₀^2
	
		rays = Ray[]
		for θ₀ in angles
			ξ₀ = cos(θ₀) / c₀
			ζ₀ = sin(θ₀) / c₀
	
			u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀Re, p₀Im, q₀Re, q₀Im]
	
			s_span = (0.0, 100e3)
		
			prob = ODEProblem(tracer!, u₀, s_span)
		
			sol = solve(prob, Tsit5(), callback = cb)
	
			r(s) = sol(s, idxs = 1)
			z(s) = sol(s, idxs = 2)
			s_max = sol.t[end]
	
			push!(rays, Ray(θ₀, r, z, s_max))
		end
		new(scn, rays)
	end
end

export Trace

# User Recipe
@userplot RayTracePlot
@recipe function plot(rtp::RayTracePlot)
	# Parse Inputs
	trc = rtp.args[1]

	# Plot Design.
	legend --> :none

	z_rng_srf = calc_bnd_range(trc.scn, :srf)
	z_rng_btm = calc_bnd_range(trc.scn, :btm)
	ex = (srf = z_rng_srf.lo, btm = z_rng_btm.hi)
	ylims --> (ex.srf, ex.btm)

	yflip := true

	# Plot Labels
	plot_title := trc.scn.name
	xguide := "Range [m]"
	yguide := "Depth [m]"

	# Ray Trace
	for ray in trc.rays
		s = range(0.0, ray.s_max, 501)
		r = ray.r.(s)
		z = ray.z.(s)
		@series r, z
	end

	# Boundaries
	for boundary in (:srf, :btm)
		bnd = getproperty(trc.scn.env, boundary)
		x = range(0.0, trc.scn.ent.rcv.x)
		z = bnd.z.(x)
		@series begin
			linecolor := :gray
			linealpha := 0.5
			fillalpha := 0.5
			fillrange := zeros(size(z)) .+ ex[boundary]
			fillstyle := :/
			fillcolor := :gray
			x, z
		end
		# @series begin
		# 	seriestype := :ribbon
		# 	[x_rng.lo, x_rng.hi], ex[boundary] .+ [0, 100]
		# end
	end
	# boundary = :btm
	# bnd = getproperty(trc.scn.env, boundary)
	# x = range(0.0, trc.scn.ent.rcv.x)
	# z = bnd.z.(x)
	# @series begin
	# 	linecolor := :black
	# 	fillrange := zeros(size(z)) .+ ex[boundary]
	# 	# fillstyle := 
	# 	x, z
	# end
	# @series begin
	# 	seriestype := :ribbon
	# 	[x_rng.lo, x_rng.hi], ex[boundary] .+ [0, 100]
	# end
end

# export RayTracePlot
export raytraceplot

# # Type Recipes
# @recipe function plot(::Type{Trace}, trc::Trace)

# end

# #Plot Recipes
# @recipe function plot(::Type{Val{:myplotrecipename}})

# end

# # Series Recipes
# @recipe function plot(::Type{Val{:myseriesrecipename}}, x, y, z)

# end
