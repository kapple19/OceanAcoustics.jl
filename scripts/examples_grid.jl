##
using OceanAcoustics
grids = example_grid.(OAC_EXAMPLE_NAMES)
p = plot_oac.(grids)
save_oac_plot.(p, :grids, OAC_EXAMPLE_NAMES)

##
using OceanAcoustics
for name ∈ OAC_EXAMPLE_NAMES[5:end]
	grid = example_grid(name)
	p = plot_oac(grid)
	display(p)
	save_oac_plot(p)
end

##
# 1 => 56
# 2 => 36
# 3 => 26, 38, 41, 36, 31
# 4 => 34, 33
# 5 => 23, 29, 22
# 6 => 20, 24, 20, 20, 23, 27, 23
# 7 => 21, 19, 19, 27, 28, 21, 19
# 8 => 21, 19, 19, 20, 23, 28, 
# 9 => 29, 23, 25
# 10 => 22, 27, 29
using OceanAcoustics
grid = example_grid(:convergence)
p = plot_oac(grid)
