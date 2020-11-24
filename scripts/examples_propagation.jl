##
using OceanAcoustics
props = example_propagation.(OAC_EXAMPLE_NAMES)
p = plot_oac.(props)
save_oac_plot.(p, :propagation, OAC_EXAMPLE_NAMES)

##
using OceanAcoustics
names_subset = OAC_EXAMPLE_NAMES[1:7]
for (n, name) ∈ enumerate(names_subset)
	print("Propagation ", n, "/", length(names_subset))
	println(": ", name)
	prop = example_propagation(name)
	p = plot_oac(prop)
	display(p)
	save_oac_plot(p, :propagation, name)
end

##
using OceanAcoustics
prop = example_propagation(:n2linear)
p = plot_oac(prop)

##
using OceanAcoustics
fld = example_field(:convergence)
prop = @time Propagation(fld)
p = plot_oac(prop)
