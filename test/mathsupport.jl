@testset "Function Interpolation" begin
	x = [0.0, 2e2, 2e3, 5e3]
	y = [1480.0, 1510.0, 1500.0, 1520.0]

	fcn = OceanAcoustics.interpolated_function(x, y)
	@test fcn(-1.0) == y[1]
	for n = 1:length(x)
		@test fcn(x[n]) == y[n]
	end
	@test fcn(x[end] + 1) == y[end]

	y = x
	x = [0.0, 5e2, 5e3, 10e3, 20e3]
	z = 1500.0 .- x/20e3*100 .+ y'/5e3*200

	fcn = OceanAcoustics.interpolated_function(x, y, z)
	for x₀ ∈ [-1.0, 0.0], y₀ ∈ [-1.0, 0.0]
		@test fcn(x₀, y₀) == z[1]
	end
	for x₀ ∈ x[end] .+ [0.0, 1.0], y₀ ∈ y[1] .- [0.0, 1.0]
		@test fcn(x₀, y₀) == z[end, 1]
	end
	for x₀ ∈ x[1] .- [0.0, 1.0], y₀ ∈ y[end] .+ [0.0, 1.0]
		@test fcn(x₀, y₀) == z[1, end]
	end
	for x₀ ∈ x[end] .+ [0.0, 1.0], y₀ ∈ y[end] .+ [0.0, 1.0]
		@test fcn(x₀, y₀) == z[end, end]
	end
	for (nx, x₀) ∈ enumerate(x), (ny, y₀) ∈ enumerate(y)
		@test fcn(x₀, y₀) == z[nx, ny]
	end
end

v₀ = 0.0:0.5:2
v₊ = v₀[2:end]
v = [-v₊; 0; v₊]
fcn_tests = [
	(
		v, v,
		(x, y) -> x^2 * y^3, # f
		(x, y) -> 2x * y^3, # ∂f_∂x
		(x, y) -> x^2 * 3y^2, # ∂f_∂y
		(x, y) -> 2y^3, # ∂²f_∂x²
		(x, y) -> 2x * 3y^2, # ∂²f_∂x∂y
		(x, y) -> 2x * 3y^2, # ∂²f_∂y∂x
		(x, y) -> x^2 * 6y # ∂²f_∂y²
	), (
		v, v,
		(x, y) -> x*exp(2y), # f
		(x, y) -> exp(2y), # ∂f_∂x
		(x, y) -> 2x*exp(2y), # ∂f_∂y
		(x, y) -> 0, # ∂²f_∂x²
		(x, y) -> 2exp(2y), # ∂²f_∂x∂y
		(x, y) -> 2exp(2y), # ∂²f_∂y∂x
		(x, y) -> 4x*exp(2y) # ∂²f_∂y²
	), (
		v₊, v₊,
		(x, y) -> y^1.5*√x, # f
		(x, y) -> y^1.5/(2*√x), # ∂f_∂x
		(x, y) -> 1.5*√(x*y), # ∂f_∂y
		(x, y) -> -y^1.5/4*x^(-3/2), # ∂²f_∂x²
		(x, y) -> 0.75*√(y/x), # ∂²f_∂x∂y
		(x, y) -> 0.75*√(y/x), # ∂²f_∂y∂x
		(x, y) -> 0.75*√(x/y) # ∂²f_∂y²
	), (
		v, v,
		(x, y) -> x + x^2/2 + 3y^3,# f
		(x, y) -> 1 + x, # ∂f_∂x
		(x, y) -> 9y^2, # ∂f_∂y
		(x, y) -> 1, # ∂²f_∂x²
		(x, y) -> 0, # ∂²f_∂x∂y
		(x, y) -> 0, # ∂²f_∂y∂x
		(x, y) -> 18y # ∂²f_∂y²
	)
]

@testset "Bivariate Partial Derivatives" for (
	xs, ys,
	f, analytic_∂f_∂x, analytic_∂f_∂y,
	analytic_∂²f_∂x²,
	analytic_∂²f_∂x∂y, analytic_∂²f_∂y∂x,
	analytic_∂²f_∂y²
	) ∈ fcn_tests

	∂f_∂x, ∂f_∂y, ∂²f_∂x², ∂²f_∂x∂y, ∂²f_∂y∂x, ∂²f_∂y² = OceanAcoustics.bivariate_derivatives(f)
	for x ∈ xs, y ∈ ys
		@test ∂f_∂x(x, y) ≈ analytic_∂f_∂x(x, y)
		@test ∂f_∂y(x, y) ≈ analytic_∂f_∂y(x, y)
		@test ∂²f_∂x²(x, y) ≈ analytic_∂²f_∂x²(x, y)
		@test ∂²f_∂x∂y(x, y) ≈ ∂²f_∂y∂x(x, y)
		@test ∂²f_∂x∂y(x, y) ≈ analytic_∂²f_∂x∂y(x, y)
		@test ∂²f_∂y∂x(x, y) ≈ analytic_∂²f_∂y∂x(x, y)
	end
end

@testset "Univariate Interpolation" begin
	
end

@testset "Bivariate Interpolation" begin
	xVec = [0, 1, 2, 5.]
	yVecVecs = [
		[0, 1, 2, 4, 6.],
		[0, 5, 15, 50.],
		[0, 100.],
		[0, 10, 20, 40.]
	]
	zVecVecs = [
		[10, 9, 10, 11, 12.],
		[5, 4, 3, 4.],
		[20, 25.],
		[1, 2, 3, 10.]
	]

	f = bivariate_interpolation(xVec, yVecVecs, zVecVecs)

	@testset "Ranged Depth Vectors" for nx = eachindex(xVec), ny = eachindex(yVecVecs[nx])
		@test f(xVec[nx], yVecVecs[nx][ny]) == zVecVecs[nx][ny]
	end

	@testset "Outside" begin
		@test zVecVecs[1][1] == f(xVec[1] - 1, yVecVecs[1][1])
		@test zVecVecs[1][1] == f(xVec[1], yVecVecs[1][1] - 1)
		@test zVecVecs[1][1] == f(xVec[1] - 1, yVecVecs[1][1] - 1)
	end

	xVec = [0, 10, 20, 50]
	zVecFcns = [
		y -> sin(y),
		y -> exp(-y),
		y -> y - 1/y,
		y -> y^2 - y^3
	]
	f = bivariate_interpolation(xVec, zVecFcns)
	@testset "Ranged Depth Functions" for (nx, x) ∈ enumerate(xVec), y ∈ LinRange(-1, 10, 5)
		@test f(x, y) == zVecFcns[nx](y)
	end

end