@testset "Examples on Ray Method" for scenario in propertynames(examples)
	scn = getproperty(examples, scenario)
	@info "Testing $(scn.name)"

	trc = Trace(scn, π/6 * range(-1.0, 1.0, 11))
	@test trc isa Trace

	p = plot()
	for ray in trc.rays
		s = range(0.0, ray.s_max, 5)
		r = ray.r.(s)
		z = ray.z.(s)
		@test all(diff(r) .> 0)
		plot!(r, z)
	end
	savefig(p, joinpath("img", "trace_" * string(scenario) * ".png"))
end