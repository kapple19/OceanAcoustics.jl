using OceanAcoustics
flds = example_field.(OAC_EXAMPLE_NAMES)
ps = plot_oac.(flds)
save_oac_plot.(ps, :fields, OAC_EXAMPLE_NAMES)
