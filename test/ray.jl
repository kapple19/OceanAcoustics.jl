@testset "Examples on Ray Method" for (scenario, scn) in pairs(examples)
	@info "Ray Method: $(scn.name)"

	trc = Trace(scn, 21)
	@test trc isa Trace

	rtp = raytraceplot(trc)
	savefig(rtp, joinpath("img", "trace_" * string(scenario) * ".png"))

	trc = Trace(scn, 101)
	r, z, TL = Field(trc)
	@test TL isa Array
	fig = heatmap(r, z, TL',
		c = cgrad(:jet, rev = true),
		title = scn.name,
		yflip = true
	)
	savefig(fig, joinpath("img", "raymethod_" * string(scenario) * ".png"))
end

@testset "Literature Replications" begin
	@info "Jensen Fig 3.16"
	@test begin
		scn = examples.n2_linear_profile
		trc = Trace(scn, [-π/4])
		r, z, TL = Field(trc)
		TL = max.(40, TL)
		TL = min.(90, TL)

		fig = heatmap(r, z, TL',
			c = cgrad(:jet, rev = true),
			title = scn.name,
			yflip = true
		)
		savefig(fig, joinpath("img", "jensenetal2011_fig_3_16.png"))

		TL isa Array
	end
end