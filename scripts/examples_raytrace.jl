## All
using OceanAcoustics
trcs = example_trace.(OAC_EXAMPLE_NAMES)
ps = plot_oac.([trc.scn.env for trc ∈ trcs])
plot_oac!.(trcs)
save_oac_plot.(ps, :raytraces, OAC_EXAMPLE_NAMES)

## One
using OceanAcoustics
name = :stacked
scn = example_scenario(name)
trc = Trace(scn)
p = plot_oac(scn.env)
plot_oac!(trc)
