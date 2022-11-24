r = unique(sort(10 * rand(20)))
z = unique(sort(10 * rand(length(r))))
alt = OceanAcoustics.Altimetry(r, z)
@testset "altimetry" begin
	@test alt.min == minimum(z)
	@test alt.max == maximum(z)
	@test all(alt.z.(r) == z)
end