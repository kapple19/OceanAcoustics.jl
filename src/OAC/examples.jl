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