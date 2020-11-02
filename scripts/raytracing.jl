## One
using OceanAcoustics
# name = OAC_EXAMPLE_NAMES[2]
name = :convergence
trc = example_trace(name)
fig = plot_oac(trc)
save_oac_plot(fig, :raytrace, name)

## All
using OceanAcoustics
trcs = example_trace.(OAC_EXAMPLE_NAMES)
figs = plot_oac.(trcs)
save_oac_plot.(figs, :raytrace, OAC_EXAMPLE_NAMES)
