##
using OceanAcoustics

##
include("../scripts/scenarios.jl")

θ₀, src, ocn, bty, ati, title = simple()

fld = Field(θ₀, src, ocn, bty, ati)

rng = range(0, ocn.R, length = 31)
dpt = range(0, ocn.Z, length = 15)
rays = Vector{Ray}(undef, 0)
for nRay = 1:length(fld.beams)
	push!(rays, fld.beams[nRay].ray)
end

## Labels
acoustic_plot()

## Limits
@show ocn.R
@show ocn.Z
acoustic_plot(ocn.R, ocn.Z)

## Boundary
rng = range(0, fld.ocn.R, length = 101)
acoustic_plot(rng, ati, ocn.Z)
acoustic_plot!(rng, bty)

## Rays
rays = Vector{Ray}(undef, 0)
for nRay = 1:length(fld.beams)
	push!(rays, fld.beams[nRay].ray)
end
acoustic_plot(rays)

## Field
rng = range(0, ocn.R, length = 31)
dpt = range(0, ocn.Z, length = 15)

acoustic_plot(rng, dpt, fld)

## Altogether
p = acoustic_plot(rng, dpt, fld)
acoustic_plot!(rng, ati)
acoustic_plot!(rng, bty)
acoustic_plot!(rays)
display(p)
