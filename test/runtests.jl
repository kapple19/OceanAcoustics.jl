using Test
using OceanAcoustics
using Plots

@testset "OceanAcoustics.jl" begin
    include("scenarios.jl")
    include("ray.jl")
end
