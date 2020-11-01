export plot_oac
export plot_oac!
export save_oac_plot

function plot_oac()
	f = Figure()
	colorscheme!(f, "light")
	hold!(f, true)
	return f
end

function plot_oac!()
	f = gcf()

	yflip!(f, true)
	xlabel!(f, "Range (m)")
	ylabel!(f, "Depth (m)")
end

function plot_oac(oac::OceanAcoustic)
	f = plot_oac()
	plot_oac!(oac)
	plot_oac!()
	return f
end

function plot_oac!(ttl::AbstractString)
	f = gcf()
	title!(f, ttl)
end

function plot_oac!(env::Environment)
	f = gcf()

	r = LinRange(
		env.Ωr.lo,
		env.Ωr.hi,
		1001
	)
	plot!(
		f, r, env.ati.z.(r),
		linecolor = 'k'
	)
	plot!(
		f, r, env.bty.z.(r),
		linecolor = color(0, 0, 0)
	)
end

# function plot_oac(env::Environment)
# 	f = plot_oac()
# 	plot_oac!(env)
# 	plot_oac!()
# 	return f
# end

function plot_oac!(ray::Ray)
	f = gcf()

	s = LinRange(ray.Ωₛ.lo, ray.Ωₛ.hi, 1001)
	plot!(
		f, ray.r.(s), ray.z.(s),
		linecolor = color(0, 0.4, 0.7)
	)
end

# function plot_oac(ray::Ray)
# 	f = plot_oac()
# 	plot_oac!(ray)
# 	plot_oac!()
# 	return f
# end

function plot_oac!(trc::Trace)
	f = gcf()
	plot_oac!(trc.scn.env)
	for nRay ∈ 1:length(trc.rays)
		plot_oac!(trc.rays[nRay])
	end
end

# function plot_oac(trc::Trace)
# 	f = plot_oac()
# 	plot_oac!(trc)
# 	plot_oac!()
# 	return f
# end

function save_oac_plot(
	f::Figure,
	filedirname::AbstractString...)
	gcf(f)
	plotsdir(filedirname...) * ".png" |> savefig
end
