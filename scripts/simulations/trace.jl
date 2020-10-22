## All
using OceanAcoustics

trcs = trace_example.(OAC_EXAMPLE_NAMES)

fs = oac_plot.(trcs)

save_oac_plot.(fs, "examples", "trace", String.(OAC_EXAMPLE_NAMES))

## One
using OceanAcoustics

example = :channel

trc = trace_example(example)

f = oac_plot(trc)

save_oac_plot(f, "examples", "trace", String(example))
