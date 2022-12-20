@testset "Scenario" begin
	for scenario in examples
		scn = getproperty(Examples, scenario)
		sp = scenarioplot(scn)
		savefig(sp, joinpath("img", "scenario", "scenario_" * string(scenario) * ".png"))
	end
end