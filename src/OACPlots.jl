export oac_plot
export save_oac_plot

function oac_plot()
	f = Figure()
	hold!(f, true)
	xlabel!(f, "Range (m)")
	ylabel!(f, "Depth (m)")
	colorscheme!(f, "light")
	return f
end

function oac_plot!()
	f = gcf()
	yflip!(f, true)
end

function oac_plot!(ttl::AbstractString)
	f = gcf()
	title!(f, ttl)
end

function oac_plot(env::Environment)
	f = oac_plot()

	r = LinRange(0.0, env.Ωr.hi, 1001)
	plot!(f, r, env.bty.z,
		linecolor = color(0, 0, 0))
	plot!(f, r, env.ati.z,
		linecolor = color(0, 0, 0))
	xlim((env.Ωr.lo, env.Ωr.hi))
	ylim((env.Ωz.lo, env.Ωz.hi))

	oac_plot!()
	return f
end

function oac_plot(trc::Trace)
	f = oac_plot(trc.sno.env)
	for ray = trc.rays
		s = LinRange(0.0, ray.S, 1001)
		plot!(f, ray.r.(s), ray.z.(s))
	end
	oac_plot!("Ray Trace: " * trc.sno.name)
	oac_plot!()
	return f
end

function save_oac_plot(f::Figure, filedirname::AbstractString...)
	gcf(f)
	savefig(plotsdir(filedirname...) * ".png")
end