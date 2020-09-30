##
using Plots
using ForwardDiff
using OceanAcoustics
using DrWatson
using Roots

##
include("../scripts/scenarios.jl")

##
dict_rays_n2linear = scenario_n2linear()
@unpack z₀, f, c, R, zBty, zAti, title = dict_rays_n2linear[:pars]

src = Source(Position(0, z₀), Signal(f))
ocn = Medium(c, R)
bty = Boundary(zBty)
ati = Boundary(zAti)

θ₀ = -acos(c(0, z₀)/c(0, 0))

beam = Beam(θ₀, src, ocn, bty, ati);

TL(s, n) = min(100, -20log10(abs(beam.b(s, n))))

##
r = 5e2
z = 5e2
s = range(0, beam.ray.S, length = 101)

##
plot(beam.ray.sol, vars = (1, 2), yaxis = :flip)
scatter!([r], [z])

##
Q(s) = (beam.ray.r(s) - r)^2 + (beam.ray.z(s) - z)^2

plot(s, Q)

##
dQ(s) = ForwardDiff.derivative(Q, s)

##
plot(s, dQ)

##
sMin = find_zero(dQ, (0, 1000))

##
plot(s, Q)
scatter!([sMin], [Q(sMin)])

nMin = sqrt(Q(sMin))
##
