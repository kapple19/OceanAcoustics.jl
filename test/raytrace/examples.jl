@testset "Examples on Ray Trace Method" for scenario in examples
	img_dir = joinpath("img", "raytrace", "examples")

	scn = getproperty(Examples, scenario)
	@info "Ray Trace Method: $(scn.name)"

	trc = RayTrace.Field(scn,
		21,
		save_field = false,
		save_trace = true
	)
	@test trc isa Trace

	rtp = raytraceplot(trc)
	scenarioplot!(scn)
	savefig(rtp,
		joinpath(img_dir, "trace_" * string(scenario) * ".png")
	)

	fld = RayTrace.Field(scn, 101)
	@test fld.r isa AbstractVector{<:AbstractFloat}
	@test fld.z isa AbstractVector{<:AbstractFloat}
	@test fld.PL isa AbstractMatrix{<:AbstractFloat}

	fig = propagationplot(fld)
	scenarioplot!(scn)
	savefig(fig,
		joinpath(img_dir, "proploss_" * string(scenario) * ".png")
	)
end