using OceanAcoustics

scn = ExampleScenarios.wavy()

trc = Trace(scn)

f = plot_oac(trc)
save_oac_plot(f, "trace_eg4")
