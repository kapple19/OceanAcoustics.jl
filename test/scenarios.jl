
# @testset "Depth" begin
# 	a, b, c, d = rand(4)
# 	fcn(x) = a * x^2 - b * exp(-x) + c * sin(x)
# 	dom = 0.0 .. 1 + d
# 	rng = fcn(dom)
# 	depth_wfcn = Depth(fcn, dom)
	
# 	x = unique(sort(dom.hi * rand(20)))
# 	z = fcn.(x)
# 	depth = Depth(x, z)

# 	@testset "Inputs" begin
# 		@testset "Scalar" begin
# 			zF64 = 10 * rand()
# 			altF64 = Depth(zF64)
# 			@test altF64(0) == zF64

# 			zI64 = rand(Int)
# 			altI64 = Depth(zI64)
# 			@test altI64(0) == Float64(zI64)
# 		end

# 		@testset "Vectors" begin
# 			@test depth.min == minimum(z)
# 			@test depth.max == maximum(z)
# 			@test depth.(x) == z
# 		end

# 		@testset "Functions" begin
# 			@test depth_wfcn.min ≈ rng.lo
# 			@test depth_wfcn.max ≈ rng.hi
# 			@test depth_wfcn.(x) ≈ z
# 		end
# 	end

# 	@testset "Exceptions" begin
# 		@test_throws DimensionMismatch Depth(x, [z; 10 * rand(3)])
# 		@test_throws NotSorted Depth([3, 2, 1], [1, 2, 3])
# 		@test_throws NotAllUnique Depth([1, 1, 2], [1, 2, 3])
# 	end
# end

@testset "Scenario" begin
	for scenario in examples
		scn = getproperty(Examples, scenario)
		sp = scenarioplot(scn)
		savefig(sp, joinpath("img", "scenario", "scenario_" * string(scenario) * ".png"))
	end
end