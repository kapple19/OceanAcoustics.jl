##
using OceanAcoustics

for scenario in propertynames(examples)[2:end]
	println("="^length(scenario |> string))
	println(scenario)

	scn = getproperty(examples, scenario)
	trc = Trace(scn, π/6 * range(-1.0, 1.0, 11))

	for ray in trc.rays
		s = range(0.0, ray.s_max, 5)
		println()
		@show ray.r.(s)
		println()
		@show ray.z.(s)
	end
end
