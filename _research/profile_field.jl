

using ProfileView

include("../src/OceanAcoustics.jl")
include("../scripts/scenarios.jl")

R = 10e3
Z = 1e3
r₀ = 0.0
z₀ = Z
f = 2e3
c₀ = 1550.
c(r, z) = c₀/sqrt(1 + 2.4z/c₀)

src = OceanAcoustics.Source(OceanAcoustics.Position(r₀, z₀), OceanAcoustics.Signal(f))
ocn = OceanAcoustics.Medium(c, R, Z)
bty = OceanAcoustics.Boundary(Z, R)
ati = OceanAcoustics.Boundary(0.0, R)

θ₀_crit = -acos(ocn.c(r₀, z₀)/ocn.c(r₀, 150.))
θ₀ = θ₀_crit*(0.1:0.05:1.1)

fld = OceanAcoustics.Field(θ₀, src, ocn, bty, ati)

#= In Julia REPL:
rng = LinRange(0.0, ocn.R, 31)
dpt = LinRange(0.0, ocn.Z, 15)
@profview fld.TL.(rng', dpt)
=#