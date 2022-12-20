@testset "Literature Replications" begin
	@info "Jensen Fig 3.16"
	@test begin
		scn = getproperty(Examples, :n2_linear_profile)
		fld = RayTrace.Field(scn, [-π/4])
		fld.PL = max.(40, fld.PL)
		fld.PL = min.(90, fld.PL)

		fig = propagationplot(fld)
		scenarioplot!(scn)
		savefig(fig, joinpath("img", "literature", "jensenetal2011_fig_3_16.png"))

		fld isa Field
	end
end