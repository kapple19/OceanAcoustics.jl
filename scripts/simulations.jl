function run_sims(sim_fcn::Function, scenarios::AbstractVector)
	for scen ∈ scenarios
		p = sim_fcn(scen)
		filename = String(Symbol(scen))
		dirname = String(Symbol(sim_fcn))
		wsave(plotsdir(dirname) * "/" * filename * ".png", p)
	end
end