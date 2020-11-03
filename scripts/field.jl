## One
using OceanAcoustics
name = OAC_EXAMPLE_NAMES[2]
scn = example_scenario(name)
fld = Field(scn)

##
rng = LinRange(scn.env.Ωr.lo, scn.env.Ωr.hi, 15)
dpt = LinRange(scn.env.Ωz.lo, scn.env.Ωz.hi, 9)
# p = @time [fld.p(r, z) for z ∈ dpt, r ∈ rng]
# TL = min.(100, -20log10.(abs.(p)))
TL = @time [fld.TL(r, z) for z ∈ dpt, r ∈ rng]

##
using GRUtils

f = Figure()
heatmap!(f, rng, dpt, TL,
	yaxis = :flip
)
yflip!(f, true)
colormap!(f, "jet")
display(f)
