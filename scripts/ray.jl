##
using OceanAcoustics

for scenario in propertynames(examples)
	scn = getproperty(examples, scenario)
	println("="^(4*length(scenario |> string)))
	println(scn.name)

	trc = RayMethodField(scn,
	11,
	save_field = false,
	save_trace = true
	)
	@show trc
end

## Just One
using OceanAcoustics

scenario = :n2_linear_profile
scn = getproperty(examples, scenario)
trc = RayMethodField(scn,
	11,
	save_field = false,
	save_trace = true
)
for ray in trc.rays
	@show ray.s_max
end

##
using OceanAcoustics

scenario = :n2_linear_profile
scn = getproperty(examples, scenario)
fld = RayMethodField(scn,
	[-π/4]
)
