## All
using OceanAcoustics
trcs = example_trace.(OAC_EXAMPLE_NAMES)
ps = plot_oac.([trc.scn.env for trc ∈ trcs])
plot_oac!.(trcs)
save_oac_plot.(ps, :raytraces, OAC_EXAMPLE_NAMES)

## Loop
# for name ∈ OAC_EXAMPLE_NAMES
ps = []
for name ∈ [:simple_convergence, :channel]
	trc = example_trace(name)
	push!(ps, plot_oac(trc.scn.env))
	plot_oac!(trc)
end
display.(ps)

## One
using OceanAcoustics
name = :channel
scn = example_scenario(name)
trc = Trace(scn)
p = plot_oac(scn.env)
plot_oac!(trc)
