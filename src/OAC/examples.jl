examples = let
	example_names = Symbol[]
	example_scenarios = Scenario[]

	env_north_atlantic = let
		ocn = Ocean(
			[0, 300, 1200, 2e3, 5e3],
			[1522, 1501, 1514, 1496, 1545.0]
		)

		Environment(ocn, 5e3)
	end

	scn_north_atlantic_convergence_zones = Scenario(
		env_north_atlantic,
		((200, 0), 70e3),
		"North Atlantic Convergence Zones"
	)

	push!(example_names, :north_atlantic_convergence_zones)
	push!(example_scenarios, scn_north_atlantic_convergence_zones)

	munk_profile = let
		f = 5e2
		z_src = 1e3
		r_rcv = 20e3

		z̃(z) = 2/1300*(z - 1300)
		ϵ = 7.37e-3
		c(r, z) = 1500(1 + ϵ*(z̃(z) - 1 + exp(-z̃(z))))

		ocn = Ocean(c)
		scn = Scenario((ocn, 5e3), ((f, z_src), r_rcv), "Munk Profile")
	end

	push!(example_names, :munk_profile)
	push!(example_scenarios, munk_profile)

	n2_linear_profile = let
		c₀ = 1550.0
		c(r, z) = c₀ / √(max(0.0, 1 + 2.4z / c₀))

		ocn = Ocean(c)

		scn = Scenario(
			(ocn, 1e3),
			((2e3, 1e3), 3.5e3),
			"n²-Linear Profile"
		)
	end

	push!(example_names, :n2_linear_profile)
	push!(example_scenarios, n2_linear_profile)

	parabolic_bathymetry = let
		r_rcv = 20e3
	
		c = 250.0
		b = 2.5e5
		z_bty(r) = 2e-3b * √(1 + r/c)
		
		btm = Bottom(z_bty, z_bty(0.0), z_bty(r_rcv))
	
		scn = Scenario(
			(c, btm),
			((2e2, 0.0), r_rcv),
			"Parabolic Bathymetry"
		)
	end
	
	push!(example_names, :parabolic_bathymetry)
	push!(example_scenarios, parabolic_bathymetry)

	(; zip(example_names, example_scenarios)...)
end

export examples