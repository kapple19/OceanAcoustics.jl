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

examples = (; zip(example_names, example_scenarios)...)

export examples