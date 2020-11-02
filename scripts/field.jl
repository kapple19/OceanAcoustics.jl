## One
using OceanAcoustics
name = OAC_EXAMPLE_NAMES[1]
scn = example_scenario(name)
fld = Field(scn)

##
rng = LinRange(scn.env.Ωr.lo, scn.env.Ωr.hi, 41)
dpt = LinRange(scn.env.Ωz.lo, scn.env.Ωz.hi, 27)
p = @time [fld.p(r, z) for r ∈ rng, z ∈ dpt]
TL = min.(100, -20log10.(abs.(p)))

##
using Plots

heatmap(rng, dpt, TL',
	yaxis = :flip,
	seriescolor = cgrad(:jet, rev = true)
)
