##
using Plots
using ForwardDiff
using OceanAcoustics
using DrWatson

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
plot(beam.ray.sol, vars = (1, 2), yaxis = :flip)

##
s = range(0, beam.ray.S, length = 101)
n = 200range(-1, 1, length = 51)

##
r(s) = beam.ray.r(s)
z(s) = beam.ray.z(s)

dr_ds(s) = ForwardDiff.derivative(r, s)
dz_ds(s) = ForwardDiff.derivative(z, s)

r(s, n) = r(s) - dz_ds(s)*n
z(s, n) = z(s) + dr_ds(s)*n

d²r_ds²(s) = ForwardDiff.derivative(dr_ds, s)
d²z_ds²(s) = ForwardDiff.derivative(dz_ds, s)

dr_ds(s, n) = dr_ds(s) - d²z_ds²(s)*n
dr_dn(s, n) = dz_ds(s)
dz_ds(s, n) = dz_ds(s) + d²r_ds²(s)*n
dz_dn(s, n) = dr_ds(s)

𝒥(s, n) = dr_ds(s, n)*dz_dn(s, n) - dr_dn(s, n)*dz_ds(s, n)

##
contour!(s, n, 𝒥, levels = [0])

##
heatmap(s, n, TL,
	seriescolor = cgrad(:jet, rev = true))
plot!(s, beam.W)

##
scatter(r.(s, n'), z.(s, n'),
	aspect_ratio = 1,
	legend = false,
	yaxis = :flip)

##
plot(r.(s), z.(s), yaxis = :flip)
rW(s) = r(s, beam.W(s))
zW(s) = z(s, beam.W(s))
plot!(rW.(s), zW.(s))

##
using IntervalRootFinding
