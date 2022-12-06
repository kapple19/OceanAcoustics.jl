@testset "Examples on Ray Method" for (scenario, scn) in pairs(examples)
	@info "Trace: $(scn.name)"

	trc = Trace(scn)
	@test trc isa Trace

	rtp = raytraceplot(trc)
	savefig(rtp, joinpath("img", "trace_" * string(scenario) * ".png"))
end