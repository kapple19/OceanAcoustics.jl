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

	bnd = Boundary(zFcn, r -> 1500)
	@testset "Function Input" for r ∈ 3LinRange(-1, 1, 11)
		@test zFcn(r) ≈ bnd.z(r)
	end

	zFcns = [r -> exp(-r); r -> cos(2r)*sin(3r)]
	cFcns = [r -> 1500 for nz ∈ zFcns]
	bnds = Boundary.(zFcns, cFcns)
	@testset "Functions Input" for r ∈ 3LinRange(-1, 1, 11), nFcn ∈ eachindex(zFcns)
			@test zFcns[nFcn](r) ≈ bnds[nFcn].z(r)
	end
end

@testset "Celerity" begin
	R = 10e3
	Z = 2e3
	cFcn(r, z) = 1500 - 100r/R + 100z/Z
	an_∂c_∂r(r, z) = -100/R
	an_∂c_∂z(r, z) = 100/Z
	an_∂²c(r, z) = 0
	SSP = OceanAcoustics.Celerity(cFcn)
	@testset "Simple Function Input" for r ∈ LinRange(0, R, 5), z ∈ LinRange(0, Z, 5)
		@test cFcn(r, z) ≈ SSP.c(r, z)
		@test an_∂c_∂r(r, z) ≈ SSP.∂c_∂r(r, z)
		@test an_∂c_∂z(r, z) ≈ SSP.∂c_∂z(r, z)
		@test an_∂²c(r, z) ≈ SSP.∂²c_∂r²(r, z)
		@test an_∂²c(r, z) ≈ SSP.∂²c_∂r∂z(r, z)
		@test an_∂²c(r, z) ≈ SSP.∂²c_∂z∂r(r, z)
		@test an_∂²c(r, z) ≈ SSP.∂²c_∂z²(r, z)
	end

	cFcn(r, z) = 1600 - 100z*(z - Z)/(Z^2/4) + 100r/R
	an_∂c_∂r(r, z) = 100/R
	an_∂c_∂z(r, z) = -400(2z - Z)/Z^2
	an_∂²c_∂r²(r, z) = 0
	an_∂²c_∂r∂z(r, z) = 0
	an_∂²c_∂z∂r(r, z) = 0
	an_∂²c_∂z²(r, z) = -800/Z^2
	SSP = OceanAcoustics.Celerity(cFcn)
	@testset "Parabolic Function Input" for r ∈ LinRange(0, R, 5), z ∈ LinRange(0, Z, 5)
		@test cFcn(r, z) ≈ SSP.c(r, z)
		@test an_∂c_∂r(r, z) ≈ SSP.∂c_∂r(r, z)
		@test an_∂c_∂z(r, z) ≈ SSP.∂c_∂z(r, z)
		@test an_∂²c_∂r²(r, z) ≈ SSP.∂²c_∂r²(r, z)
		@test an_∂²c_∂r∂z(r, z) ≈ SSP.∂²c_∂r∂z(r, z)
		@test an_∂²c_∂z∂r(r, z) ≈ SSP.∂²c_∂z∂r(r, z)
		@test an_∂²c_∂z²(r, z) ≈ SSP.∂²c_∂z²(r, z)
	end
end
