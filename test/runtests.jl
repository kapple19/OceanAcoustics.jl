## Preamble
using OceanAcoustics
using Test

## Function Interpolation
@testset "Augmenting Mathematics Functions" begin
    include("augmentingmaths.jl")
end

## Environment
@testset "Ocean Acoustics Propagation" begin
    include("environment.jl")
    # include("scenario.jl")
    # include("propagation.jl")
end

## Plots
# @testset "Ocean Acoustic Plot" begin
#     include("plotting.jl")
# end

# ## Reflection
# @testset "BoundaryReflection" begin
#     for θ_bnd = deg2rad.(-85.0:5.0:85.0)
#         for θ_inc = deg2rad.(-85.0:5.0:85.0)
#             t_inc = [cos(θ_inc), sin(θ_inc)]
#             t_bnd = [cos(θ_bnd), sin(θ_bnd)]
#             t_rfl = OceanAcoustics.boundary_reflection(t_inc, t_bnd)

#             @test abs(rad2deg(θ_inc + θ_bnd)) == abs(rad2deg(θ_inc + θ_bnd))
#         end
#     end
# end

# ## Position
# @testset "Position" begin
#     for x ∈ range(0.0, 10e3, length = 51), y ∈ range(0.0, 5e3, length = 31)
#         pos = Position(x, y)
#         @test pos.r == x
#         @test pos.z == y
#     end
# end

# ## Signal
# @testset "Signal" begin
#     for f ∈ 10.0.^(0:0.5:5)
#         sig = Signal(f)
#         @test f == sig.f
#     end
# end

# ## Source
# @testset "Source" begin
#     for x ∈ range(0.0, 10e3, length = 21), y ∈ range(0.0, 5e3, length = 11)
#         pos = Position(x, y)
#         for f ∈ 10.0.^(0:5)
#             sig = Signal(f)
#             src = Source(pos, sig)
#             @test src.pos.r == x
#             @test src.pos.z == y
#             @test src.sig.f == f
#         end
#     end
# end

# ## Boundary

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

# ## Ray
# include("raytests.jl")

## Beam

## Receiver

## Proximity

## Field
