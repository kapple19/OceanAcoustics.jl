## One
using OceanAcoustics
name = OAC_EXAMPLE_NAMES[1]
scn = example_scenario(name)
fld = Field(scn)

##
rng = LinRange(scn.env.Ωr.lo, scn.env.Ωr.hi, 7)
dpt = LinRange(scn.env.Ωz.hi, scn.env.Ωz.hi, 5)
p = [fld.p(r, z) for r ∈ rng, z ∈ dpt]

##
using Plots

contourf(
	rng,
	dpt,
	abs.(p)
)