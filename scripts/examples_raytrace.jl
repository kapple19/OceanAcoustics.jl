## All
using OceanAcoustics
trcs = example_trace.(OAC_EXAMPLE_NAMES)
ps = plot_oac.([trc.scn.env for trc ∈ trcs])
plot_oac!.(trcs)
save_oac_plot.(ps, :raytraces, OAC_EXAMPLE_NAMES)

## One
using OceanAcoustics
name = :convergence
trc = example_trace(name)
p = plot_oac(trc)
