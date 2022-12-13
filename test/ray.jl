@testset "Examples on Ray Method" for (scenario, scn) in pairs(examples)
	@info "Ray Method: $(scn.name)"

	trc = RayMethods.Field(scn,
		21,
		save_field = false,
		save_trace = true
	)
	@test trc isa Trace

	rtp = raytraceplot(trc)
	scenarioplot!(scn)
	savefig(rtp, joinpath("img", "trace_" * string(scenario) * ".png"))

	fld = RayMethods.Field(scn, 101)
	@test fld.r isa AbstractVector{<:AbstractFloat}
	@test fld.z isa AbstractVector{<:AbstractFloat}
	@test fld.PL isa AbstractMatrix{<:AbstractFloat}

	fig = propagationplot(fld)
	scenarioplot!(scn)
	savefig(fig, joinpath("img", "raymethod_" * string(scenario) * ".png"))
end

@testset "Literature Replications" begin
	@info "Jensen Fig 3.16"
	@test begin
		scn = examples.n2_linear_profile
		fld = RayMethods.Field(scn, [-π/4])
		fld.PL = max.(40, fld.PL)
		fld.PL = min.(90, fld.PL)

		fig = propagationplot(fld)
		scenarioplot!(scn)
		savefig(fig, joinpath("img", "jensenetal2011_fig_3_16.png"))

		fld isa Field
	end
end