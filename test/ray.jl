@testset "Examples on Ray Method" for (scenario, scn) in pairs(examples)
	@info "Trace: $(scn.name)"

	trc = Trace(scn, (π/2 - π/7) * range(-1.0, 1.0, 21))
	@test trc isa Trace

	rtp = raytraceplot(trc)
	savefig(rtp, joinpath("img", "trace_" * string(scenario) * ".png"))
end