r = unique(sort(10 * rand(20)))
z = unique(sort(10 * rand(length(r))))
alt = Depth(r, z)
@testset "depth" begin
	@testset "inputs" begin
		@testset "vectors" begin
			@test alt.min == minimum(z)
			@test alt.max == maximum(z)
			@test all(alt.z.(r) == z)
		end

		@testset "scalar" begin
			zF64 = 10 * rand()
			altF64 = Depth(zF64)
			@test altF64.z(0) == zF64

			zI64 = rand(Int)
			altI64 = Depth(zI64)
			@test altI64.z(0) == Float64(zI64)
		end
	end

	@testset "exceptions" begin
		@test_throws DimensionMismatch Depth(r, [z; 10 * rand(3)])
		@test_throws NotSorted Depth([3, 2, 1], [1, 2, 3])
		@test_throws NotAllUnique Depth([1, 1, 2], [1, 2, 3])
	end
end