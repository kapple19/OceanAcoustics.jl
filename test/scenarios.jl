r = unique(sort(10 * rand(20)))
z = unique(sort(10 * rand(length(r))))
alt = Altimetry(r, z)
@testset "altimetry" begin
	@testset "inputs" begin
		@test alt.min == minimum(z)
		@test alt.max == maximum(z)
		@test all(alt.z.(r) == z)
	end

	@testset "exceptions" begin
		@test_throws DimensionMismatch Altimetry(r, [z; 10 * rand(3)])
		@test_throws NotSorted Altimetry([3, 2, 1], [1, 2, 3])
		@test_throws NotAllUnique Altimetry([1, 1, 2], [1, 2, 3])
	end
end