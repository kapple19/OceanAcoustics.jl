## Preamble
using OceanAcoustics
using Test

## Function Interpolation
@testset "FunctionInterpolation" begin
    x = [0, 2e2, 2e3, 5e3]
    y = [1480.0, 1510.0, 1500.0, 1520.0]

    fcn = OceanAcoustics.interpolated_function(x, y)
    @test fcn(-1) == y[1]
    for n = 1:length(x)
        @test fcn(x[n]) == y[n]
    end
    @test fcn(x[end] + 1) == y[end]

    y = x
    x = [0, 5e2, 5e3, 10e3, 20e3]
    z = 1500.0 .- x/20e3*100 .+ y'/5e3*200

    fcn = OceanAcoustics.interpolated_function(x, y, z)
    for x₀ ∈ [-1, 0], y₀ ∈ [-1, 0]
        @test fcn(x₀, y₀) == z[1]
    end
    for x₀ ∈ x[end] .+ [0, 1], y₀ ∈ y[1] .- [0, 1]
        @test fcn(x₀, y₀) == z[end, 1]
    end
    for x₀ ∈ x[1] .- [0, 1], y₀ ∈ y[end] .+ [0, 1]
        @test fcn(x₀, y₀) == z[1, end]
    end
    for x₀ ∈ x[end] .+ [0, 1], y₀ ∈ y[end] .+ [0, 1]
        @test fcn(x₀, y₀) == z[end, end]
    end
    for (nx, x₀) ∈ enumerate(x), (ny, y₀) ∈ enumerate(y)
        @test fcn(x₀, y₀) == z[nx, ny]
    end
end

## Reflection
@testset "BoundaryReflection" begin
    for θ_bnd = deg2rad.(-85:5:85)
        for θ_inc = deg2rad.(-85:5:85)
            t_inc = [cos(θ_inc), sin(θ_inc)]
            t_bnd = [cos(θ_bnd), sin(θ_bnd)]
            t_rfl = OceanAcoustics.boundary_reflection(t_inc, t_bnd)

            @test abs(rad2deg(θ_inc + θ_bnd)) == abs(rad2deg(θ_inc + θ_bnd))
        end
    end
end

## Position
@testset "Position" begin
    for x ∈ range(0, 10e3, length = 51), y ∈ range(0, 5e3, length = 31)
        pos = Position(x, y)
        @test pos.r == x
        @test pos.z == y
    end
end

## Signal
@testset "Signal" begin
    for f ∈ 10.0.^(0:0.5:5)
        sig = Signal(f)
        @test f == sig.f
    end
end

## Source
@testset "Source" begin
    for x ∈ range(0, 10e3, length = 21), y ∈ range(0, 5e3, length = 11)
        pos = Position(x, y)
        for f ∈ 10.0.^(0:5)
            sig = Signal(f)
            src = Source(pos, sig)
            @test src.pos.r == x
            @test src.pos.z == y
            @test src.sig.f == f
        end
    end
end

## Boundary

## Medium
@testset "Medium" begin
    R = 1e3
    Z = 1e3

    cVal = 1500
    ocn = Medium(cVal, R, Z)
    @test ocn.c(-1, 0) == cVal
    @test ocn.c(R/2, 0) == cVal
    @test ocn.c(R + 1, 0) == cVal

    zVec = range(0, Z, length = 3)
    cVec = [1550, 1500, 1600]
    ocn = Medium(zVec, cVec, R, Z)
    @test ocn.c(0, -1) == cVec[1]
    @test ocn.c(0, zVec[1]) == cVec[1]
    @test ocn.c(0, zVec[2]) == cVec[2]
    @test ocn.c(0, zVec[3]) == cVec[3]
    @test ocn.c(0, Z + 1) == cVec[3]
end

## Propagation

## Ray

## Beam

## Receiver

## Proximity

## Field
