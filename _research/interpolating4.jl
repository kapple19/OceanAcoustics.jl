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

beams = Beam.(θ₀, src, ocn, bty, ati);

nothing

##
plot(beams.ray.sol, vars = (1, 2),
	yaxis = :flip)

##
dr_ds(s) = ForwardDiff.derivative(s -> beam.ray.r(s), s)
dz_ds(s) = ForwardDiff.derivative(s -> beam.ray.z(s), s)
rng(s, n) = beam.ray.r(s) - dz_ds(s)*n
dpt(s, n) = beam.ray.z(s) + dr_ds(s)*n

##
function closest_points(r, z, beam)
	Q(s) = (beam.ray.r(s) - r)^2 + (beam.ray.z(s) - z)^2
	dQ(s) = ForwardDiff.derivative(Q, s)
	sMins = find_zeros(dQ, 0, beam.ray.S)
	d²Q(s) = ForwardDiff.derivative(dQ, s)
	min_cond(s) = d²Q(s) > 0 && beam.W(s) > sqrt(Q(s))
	min_cond.(sMins)
	filter!(min_cond, sMins)
	return sMins, sqrt.(Q.(sMins))
end

##
function add_to_field!(p, nr, r, nz, z, beam)
	sMins, nMins = closest_points(r, z, beam)
	for i = 1:length(sMins)
		p[nr, nz] += beam.b(sMins[i], nMins[i])
	end
end

##
r = 3.2e3
z = 9.5e2

##
pair_ = closest_points(r, z, beams)

##
Z = 1e3
ranges = range(0, R, length = 101)
depths = range(0, Z, length = 51)
p = zeros(Complex, length(ranges), length(depths))
# for beam = beams
	for (nr, r) ∈ enumerate(ranges), (nz, z) ∈ enumerate(depths)
		add_to_field!(p, nr, r, nz, z, beams)
	end
# end

##
heatmap(ranges, depths, min.(100., -20log10.(abs.(p'))),
yaxis = :flip,
seriescolor = cgrad(:jet, rev = true))
