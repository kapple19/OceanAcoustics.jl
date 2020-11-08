using OceanAcoustics
grids = example_grid.(OAC_EXAMPLE_NAMES)
p = plot_oac.(grids)
save_oac_plot.(p, :grids, OAC_EXAMPLE_NAMES)

##