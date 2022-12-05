using Test
using OceanAcoustics
using Plots

@testset "OceanAcoustics.jl" begin
    include("scenarios.jl")
    include("ray.jl")
end

@info "Remember to inspect output plots in `./test/img/`."