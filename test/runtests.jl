# Preamble
using OceanAcoustics
using Test

## Augmenting Mathematics
@testset "Augmenting Mathematics Functions" begin
	include("mathsupport.jl")
end

## Environment
@testset "Ocean Acoustics Propagation" begin
	include("environment.jl")
end
