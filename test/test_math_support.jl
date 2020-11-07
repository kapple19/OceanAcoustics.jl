@testset "interpolated_function" begin
	x = [0, 10]
	y = [100, 200]
	f = OceanAcoustics.interpolated_function(x, y)
	@test f.(x) == y
end