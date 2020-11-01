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

function plot_oac!(ray::Ray)
	f = gcf()

	s = LinRange(ray.Ωₛ.lo, ray.Ωₛ.hi, 1001)
	plot!(
		f, ray.r.(s), ray.z.(s),
		linecolor = color(0, 0.4, 0.7)
	)
end

function plot_oac!(trc::Trace)
	f = gcf()
	plot_oac!(trc.scn.env)
	plot_oac!(trc.scn.name)
	for nRay ∈ 1:length(trc.rays)
		plot_oac!(trc.rays[nRay])
	end
end

function save_oac_plot(
	f::Figure,
	file_dir_name::AbstractString...)
	gcf(f)
	plotsdir(file_dir_name...) * ".png" |> savefig
end

function save_oac_plot(
	f::Figure,
	file_dir_name...)
	fdn = [String(dir) for dir ∈ file_dir_name]
	save_oac_plot(f, fdn...)
end