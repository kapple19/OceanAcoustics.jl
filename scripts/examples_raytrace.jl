using OceanAcoustics
trcs = example_trace.(OAC_EXAMPLE_NAMES)
ps = plot_oac.(trcs)
save_oac_plot.(ps, :raytrace, OAC_EXAMPLE_NAMES)
