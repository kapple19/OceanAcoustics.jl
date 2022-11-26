##
using OceanAcoustics

for scenario in propertynames(examples)[2:end]
	println()
	println(scenario)

	scn = getproperty(examples, scenario)
	rays = Trace(scn)

	for ray in rays
		s = range(0.0, ray.s_max, 5)
		println()
		@show ray.r.(s)
		println()
		@show ray.z.(s)
	end
end