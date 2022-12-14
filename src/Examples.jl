module Examples
using ..OACBase

env_north_atlantic = let
	ocn = Ocean(
		[0, 300, 1200, 2e3, 5e3],
		[1522, 1501, 1514, 1496, 1545.0]
	)

	Environment(ocn, 5e3)
end

north_atlantic_convergence_zones = Scenario(
	env_north_atlantic,
	((200, 0), 70e3),
	"North Atlantic Convergence Zones"
)

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

parabolic_bathymetry = let
	r_rcv = 20e3

	c = 250.0
	b = 2.5e5
	z_bty(r) = 2e-3b * √(1 + r/c)
	
	btm = Bottom(z_bty)

	scn = Scenario(
		(c, btm),
		((2e2, 0.0), r_rcv),
		"Parabolic Bathymetry"
	)
end

lloyd_mirror = let
	scn = Scenario(
		(1500, 500),
		((150, 25), 500),
		"Lloyd Mirror"
	)
end

export north_atlantic_convergence_zones
export munk_profile
export n2_linear_profile
export parabolic_bathymetry
export lloyd_mirror
end # module Examples

examples = (
	north_atlantic_convergence_zones = Examples.north_atlantic_convergence_zones,
	munk_profile = Examples.munk_profile,
	n2_linear_profile = Examples.n2_linear_profile,
	parabolic_bathymetry = Examples.parabolic_bathymetry,
	lloyd_mirror = Examples.lloyd_mirror
)

export examples