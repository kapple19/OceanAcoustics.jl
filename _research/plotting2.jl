##
using OceanAcoustics

##
include("../scripts/scenarios.jl")

θ₀, src, ocn, bty, ati, title = simple()

fld = Field(θ₀, src, ocn, bty, ati)

##
rng = range(0, ocn.R, length = 31)
dpt = range(0, ocn.Z, length = 15)
rays = Vector{Ray}(undef, 0)
for nRay = 1:length(fld.beams)
	push!(rays, fld.beams[nRay].ray)
end

p = acoustic_plot(rng, dpt, fld)
acoustic_plot!(rng, ati)
acoustic_plot!(rng, bty)
acoustic_plot!(rays)
display(p)
