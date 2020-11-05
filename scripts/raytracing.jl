## One
using OceanAcoustics
# name = OAC_EXAMPLE_NAMES[2]
name = :convergence
trc = example_trace(name)
p = plot_oac(trc)

## All
using OceanAcoustics
trcs = example_trace.(OAC_EXAMPLE_NAMES)
ps = plot_oac.(trcs)
save_oac_plot.(ps, :raytrace, OAC_EXAMPLE_NAMES)

## Loop
using OceanAcoustics
for name ∈ OAC_EXAMPLE_NAMES
	println(name)
	trc = example_trace(name)
	p = plot_oac(trc)
	save_oac_plot(p, :raytrace, name)
end
