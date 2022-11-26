struct Ray
	θ₀
	r
	z
	s_max
end

struct Trace
	scn::Scenario
	rays::Vector{Ray}

	function Trace(scn::Scenario, angles::AbstractVector{<:AbstractFloat})
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
	
		function tracer!(du, u, params, s)
			# r, z, ξ, ζ, τ, pRe, pIm, qRe, qIm = u
			r, z, ξ, ζ = u
			
			c = cel(r, z)
			∂cr = ∂c_∂r(r, z)
			∂cz = ∂c_∂z(r, z)
			# ∂²cnn = ∂²c_∂n²(r, z, ξ, ζ)
	
			du[1] = dr = c * ξ
			du[2] = dz = c * ζ
			du[3] = dξ = -∂cr / c^2
			du[4] = dζ = -∂cz / c^2
			# du[5] = dτ = 1/c
			# dpRe = -∂²cnn / c^2 * qRe
			# dpIm = -∂²cnn / c^2 * qIm
			# dqRe = c * pRe
			# dqIm = c * pIm
	
			# du = [dr, dz, dξ, dζ, dτ, dpRe, dpIm, dqRe, dqIm]
			# du = [dr, dz, dξ, dζ]
		end
	
		f = scn.ent.src.f
	
		r₀ = 0.0
		z₀ = scn.ent.src.z
		c₀ = cel(r₀, z₀)
		# τ₀ = 0.0
		# p₀Re = 1.0
		# p₀Im = 0.0
		# λ₀ = c₀ / f
		# W₀ = 25λ₀
		# q₀Re = 0.0
		# q₀Im = pi * f * W₀^2
	
		rays = Ray[]
		for θ₀ in angles
			ξ₀ = cos(θ₀) / c₀
			ζ₀ = sin(θ₀) / c₀
	
			# u₀ = [r₀, z₀, ξ₀, ζ₀, τ₀, p₀Re, p₀Im, q₀Re, q₀Im]
			u₀ = [r₀, z₀, ξ₀, ζ₀]
	
			s_span = (0.0, 100e3)
		
			prob = ODEProblem(tracer!, u₀, s_span)
		
			sol = solve(prob, Tsit5())
	
			r(s) = sol(s, idxs = 1)
			z(s) = sol(s, idxs = 2)
			s_max = sol.t[end]
	
			push!(rays, Ray(θ₀, r, z, s_max))
		end
		new(scn, rays)
	end
end

export Trace