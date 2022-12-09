@testset "Examples on Ray Method" for (scenario, scn) in pairs(examples)
	@info "Ray Method: $(scn.name)"

	trc = RayMethodField(scn,
		21,
		save_field = false,
		save_trace = true
	)
	@test trc isa Trace

	rtp = raytraceplot(trc)
	scenarioplot!(scn)
	savefig(rtp, joinpath("img", "trace_" * string(scenario) * ".png"))

	fld = RayMethodField(scn, 101)
	@test fld.r isa AbstractVector{<:AbstractFloat}
	@test fld.z isa AbstractVector{<:AbstractFloat}
	@test fld.TL isa AbstractMatrix{<:AbstractFloat}
	fig = heatmap(fld.r, fld.z, fld.TL',
		c = cgrad(:jet, rev = true),
		title = scn.name,
		yflip = true
	)
	scenarioplot!(scn)
	savefig(fig, joinpath("img", "raymethod_" * string(scenario) * ".png"))
end

@testset "Literature Replications" begin
	@info "Jensen Fig 3.16"
	@test begin
		scn = examples.n2_linear_profile
		fld = RayMethodField(scn, [-π/4])
		fld.TL = max.(40, fld.TL)
		fld.TL = min.(90, fld.TL)

		fig = heatmap(fld.r, fld.z, fld.TL',
			c = cgrad(:jet, rev = true),
			title = scn.name,
			yflip = true
		)
		scenarioplot!(scn)
		savefig(fig, joinpath("img", "jensenetal2011_fig_3_16.png"))

		fld isa Field
	end
end