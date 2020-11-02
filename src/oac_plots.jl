export plot_oac
export plot_oac!
export save_oac_plot

gridpoints(Ω::Interval, N::Integer = 1001) = LinRange(Ω.lo, Ω.hi, N)

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

# function plot_oac!(bnd::Boundary)
# 	f = gcf()
# 	r = gridpoints(bnd.Ωr)
# 	plot!(f, r, bnd.z.(r), linecolor = 'k')
# end

# function plot_oac!(med::Medium)
# 	f = gcf()
# 	r, z = gridpoints.([med.Ωr, med.Ωz])
# 	colormap!(f, "bluescale")
# 	heatmap!(f, r, z, med.SSP.c.(r', z))
# end

function plot_oac!(env::Environment)
	f = gcf()

	function c(r, z)
		if env.ati.z(r) < z < env.bty.z(r)
			return env.ocn.SSP.c.(r, z)
		else
			return NaN
		end
	end

	r, z = gridpoints.([env.Ωr, env.Ωz])
	colormap!(f, "winter")
	contourf!(
		f, r, z, c.(r', z),
		levels = 8, majorlevels = 1
	)
end

function plot_oac!(ray::Ray)
	f = gcf()
	
	s = gridpoints(ray.Ωₛ)
	plot!(
		f, ray.r.(s), ray.z.(s),
		linecolor = color(0, 0, 0.5)
	)
	# color(0, 0.4, 0.7)
end

function plot_oac!(trc::Trace)
	f = gcf()
	plot_oac!(trc.scn.env)
	plot_oac!.(trc.rays)
	# for nRay ∈ 1:length(trc.rays)
	# 	plot_oac!(trc.rays[nRay])
	# end
	plot_oac!("Ray Trace: " * trc.scn.name)
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
