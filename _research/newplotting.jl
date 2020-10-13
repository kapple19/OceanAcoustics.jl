##
using OceanAcoustics
using GRUtils
# using Plots

#
include("../scripts/scenarios.jl")

#
θ₀, src, ocn, bty, ati, title = n2linear()

#
fld = Field(θ₀, src, ocn, bty, ati)

##
rng = range(0, ocn.R, length = 51)
dpt = range(0, ocn.Z, length = 25)

# TL = [fld.TL(r, z) for z ∈ dpt, r ∈ rng]
TL = fld.TL.(rng', dpt)

##
p = GRUtils.heatmap(TL)
yflip(true)
cmap = colormap("jet")

##
r = src.pos.r
z = src.pos.z .- [1e-12, 1e-13, 0.0]
@show fld.TL.(r, z)
@show fld.p.(r, z)
@show abs.(fld.p.(r, z))
@show abs.(fld.p.(r, z))/4π
@show -20log10.(abs.(fld.p.(r, z))/4π)
@show min.(100.0, -20log10.(abs.(fld.p.(r, z))/4π))

##
