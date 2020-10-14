function run_sims(sim_fcn::Function, scenarios::AbstractVector)
	for scen ∈ scenarios
		f = sim_fcn(scen)
		filename = String(Symbol(scen))
		dirname = String(Symbol(sim_fcn))
		savefig(plotsdir(dirname) * "\\" * filename * ".png", f)
	end
end