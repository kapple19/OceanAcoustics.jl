##
using OceanAcoustics

for scenario in propertynames(examples)[2:end]
	scn = getproperty(examples, scenario)
	println("="^(4*length(scenario |> string)))
	println(scn.name)

	trc = Trace(scn, π/6 * range(-1.0, 1.0, 11))

	for ray in trc.rays
		s = range(0.0, ray.s_max, 5)
		@show ray.r.(s)
		@show ray.z.(s)
	end
end
