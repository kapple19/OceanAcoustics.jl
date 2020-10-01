## Preamble
using OceanAcoustics
using Plots
using DrWatson

## Scenarios
include("../scripts/scenarios.jl")

## 
dict_rays_n2linear = scenario_n2linear()
@unpack z₀, f, c, R, zBty, zAti, title = dict_rays_n2linear[:pars]

Z = 1e3

src = Source(Position(0, z₀), Signal(f))
ocn = Medium(c, R, Z)
bty = Boundary(zBty)
ati = Boundary(zAti)

θ₀ = -acos(c(0, z₀)/c(0, 0))

##
rcv = Receiver(range(0, R, length = 101), range(0, Z, length = 51))

fld = Field([θ₀], src, rcv, ocn, bty, ati)

##
TL = min.(100, -20log10.(abs.(p.p)))
heatmap(rcv.r, rcv.z, TL', yaxis = :flip, seriescolor = cgrad(:jet, rev = true))

##