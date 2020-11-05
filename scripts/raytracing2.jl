## One
using OceanAcoustics
name = :flat
trc = example_trace(name)

##
p = plot_oac(trc.scn.env.bty)
p = plot_oac(trc.scn.env.ocn)
p = plot_oac(trc.scn.env)
p = plot_oac(trc)

##
using IntervalArithmetic
using Plots

SSP_COLORMAP = :winter
SSP_LEVEL_COLORS = :ice
gridpoints(Ω::Interval, N::Integer = 1001) = LinRange(Ω.lo, Ω.hi, N)

env = trc.scn.env

p = plot(yaxis = :flip, legend = false)
c(r, z) = env.ati.z(r) ≤ z ≤ env.bty.z(r) ? env.ocn.SSP.c(r, z) : NaN
r, z = gridpoints.([env.Ωr, env.Ωz])
heatmap!(
	r, z, c,
	seriescolor = SSP_COLORMAP,
	colorbar = true
)
for ray ∈ trc.rays
	s = gridpoints(ray.Ωs)
	plot!(ray.r.(s), ray.z.(s),
		color = :black)
end
plot!(title = "Ray Trace: " * trc.scn.name)
xlabel!("Range (m)")
ylabel!("Depth (m)")