using OceanAcoustics
name = :n2linear
fld = example_field(name)
ps = plot_oac(fld)
save_oac_plot(ps, :fields, name)

##
using OceanAcoustics
using IntervalArithmetic
using ProgressMeter
using Plots
fld = example_field(:n2linear)

gridpoints(Ω::Interval, N::Integer = 1001) = LinRange(Ω.lo, Ω.hi, N)
# r, z = gridpoints.([fld.scn.env.Ωr, fld.scn.env.Ωz], [1001, 501])
r, z = gridpoints.([fld.scn.env.Ωr, fld.scn.env.Ωz], [11, 7])

TL(r, z) = fld.scn.env.ati.z(r) ≤ z ≤ fld.scn.env.bty.z(r) ? min(100.0, fld.TL(r, z)) : NaN

TLgrid = @showprogress 1 "Transmission Loss Grid " [TL(r′, z′) for z′ ∈ z, r′ ∈ r]

heatmap(
	r, z, TLgrid,
	seriescolor = cgrad(:jet, rev = true),
	colorbar = true,
	colorbar_title = "Transmission Loss (dB)",
	title = fld.scn.name,
	yaxis = (:flip, "Depth (m)"),
	xaxis = "Range (m)"
)
