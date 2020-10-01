using OceanAcoustics
using Test

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
