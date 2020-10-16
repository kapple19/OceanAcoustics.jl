@testset "BoundaryReflection" begin
    for θ_bnd = deg2rad.(-85.0:5.0:85.0)
        for θ_inc = deg2rad.(-85.0:5.0:85.0)
            t_inc = [cos(θ_inc), sin(θ_inc)]
            t_bnd = [cos(θ_bnd), sin(θ_bnd)]
            t_rfl = OceanAcoustics.boundary_reflection(t_inc, t_bnd)

            @test abs(rad2deg(θ_inc + θ_bnd)) == abs(rad2deg(θ_inc + θ_bnd))
        end
    end
end

@testset "Boundary" begin
	zFcn(r) = r^2 - sin(r)
	dz_dr(r) = 2r - cos(r)

	bnd = Boundary(zFcn)
	@testset "Function" for r ∈ 3LinRange(-1, 1, 11)
		@test dz_dr(r) ≈ bnd.dz_dr(r)
	end

	zFcns = [r -> exp(-r); r -> cos(2r)*sin(3r)]
	dz_drs = [r -> -exp(-r); r -> 3cos(2r)*cos(3r) - 2sin(2r)*sin(3r)]
	bnds = Boundary.(zFcns)
	@testset "Function" for r ∈ 3LinRange(-1, 1, 11), nFcn ∈ eachindex(zFcns)
			@test bnds[nFcn].dz_dr(r) ≈ dz_drs[nFcn](r)
	end
end

# ## Medium
# @testset "Medium" begin
#     R = 1e3
#     Z = 1e3

#     cVal = 1500.0
#     ocn = Medium(cVal, R, Z)
#     @test ocn.c(-1, 0) == cVal
#     @test ocn.c(R/2, 0) == cVal
#     @test ocn.c(R + 1, 0) == cVal

#     zVec = range(0.0, Z, length = 3)
#     cVec = [1550.0, 1500.0, 1600.0]
#     ocn = Medium(zVec, cVec, R, Z)
#     @test ocn.c(0.0, -1.0) == cVec[1]
#     @test ocn.c(0.0, zVec[1]) == cVec[1]
#     @test ocn.c(0.0, zVec[2]) == cVec[2]
#     @test ocn.c(0.0, zVec[3]) == cVec[3]
#     @test ocn.c(0.0, Z + 1.0) == cVec[3]
# end