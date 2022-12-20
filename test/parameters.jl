@testset "Parameters" begin
	Depth = OceanAcoustics.OACBase.Depth
	@testset "Depth" for _ in 1:10
		a, b, c, d = rand(4)
		fcn(r) = a * r^2 - b * exp(-r) + c * sin(r)
		depth_wfcn = Depth(fcn)
		r = unique(sort((1 + d)*rand(20)))
		z = fcn.(r)

		@testset "Inputs" begin
			@testset "Function" begin
				@test depth_wfcn.(r) ≈ z
			end

			@testset "Scalar" begin
				z_scl = 100(rand() - 0.5)
				depth_wscl = Depth(z_scl)
				@test all(depth_wscl.(rand(10)) .== z_scl)
			end

			@testset "Vectors" begin
				depth_wvec = Depth(r, z)
				@test depth_wvec.(r) == z
			end

			@testset "Splat" begin
				depth_wvec_splat = Depth((r, z)...)
				@test depth_wvec_splat.(r) == z
			end
		end

		@testset "Exceptions" begin
			@test_throws DimensionMismatch Depth(r, [z; 10 * rand(3)])
			@test_throws NotSorted Depth([3, 2, 1], [1, 2, 3])
			@test_throws NotAllUnique Depth([1, 1, 2], [1, 2, 3])
		end
	end
end