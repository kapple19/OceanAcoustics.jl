using Plots:
Plot,
plot,
xlabel!,
ylabel!,
plot!,
heatmap!,
cgrad,
contourf!,
savefig
using DrWatson: plotsdir
using ProgressMeter: @showprogress

export plot_oac
export plot_oac!
export save_oac_plot

SSP_COLORMAP = :winter
SSP_LEVEL_COLOURS = :ice
DEFAULT_GRID_N = 1001

gridpoints(Ω_lo::Real, Ω_hi::Real, N::Integer = DEFAULT_GRID_N) = LinRange(Ω_lo, Ω_hi, N)

gridpoints(Ω::Tuple{Rlo, Rhi}, N::Integer = DEFAULT_GRID_N) where {Rlo <: Real, Rhi <: Real} = gridpoints(Ω[1], Ω[2], N)

gridpoints(Ω::Interval, N::Integer = DEFAULT_GRID_N) = gridpoints(Ω.lo, Ω.hi, N)

function plot_oac()
	p = plot(yaxis = :flip, legend = false)
end

function plot_oac!()
	xlabel!("Range (m)")
	ylabel!("Depth (m)")
end

function plot_oac!(bnd::Boundary, Ωr)
	r = gridpoints(Ωr)
	plot!(
		r, bnd.z,
		color = :black
	)
end

function plot_oac!(bnd::Boundary)
	plot_oac!(bnd, bnd.Ωr)
end

function plot_oac!(med::Medium, Ωr, Ωz)
	r, z = gridpoints.([Ωr, Ωz])
	if med.input_data_type ≠ :Constant
		heatmap!(
			r, z, med.SSP.c,
			seriescolor = SSP_COLORMAP,
			linewidth = 0,
			colorbar = true,
			colorbar_title = "Sound Speed Profile (m/s)"
		)
	end
end

function plot_oac!(med::Medium)
	plot_oac!(med, med.Ωr, med.Ωz)
end

function plot_oac!(env::Environment)
	c(r, z) = env.ati.z(r) ≤ z ≤ env.bty.z(r) ? env.ocn.SSP.c(r, z) : NaN
	r, z = gridpoints.([env.Ωr, env.Ωz])
	if env.ocn.input_data_type ≠ :Constant
		contourf!(
			r, z, c,
			seriescolor = SSP_COLORMAP,
			linewidth = 0,
			colorbar = true,
			colorbar_title = "Sound Speed Profile (m/s)"
		)
	end
end

function plot_oac!(ray::Ray)
	s = gridpoints(ray.Ωs)
	plot!(ray.r.(s), ray.z.(s),
		color = :navyblue)
end

function plot_oac!(trc::Trace)
	plot_oac!(trc.scn.env)
	plot_oac!.(trc.rays)
	plot!(title = "Ray Trace: " * trc.scn.name)
end

function plot_oac!(fld::Field)
	# r, z = gridpoints.([fld.scn.env.Ωr, fld.scn.env.Ωz], [1001, 801])
	r, z = gridpoints.([fld.scn.env.Ωr, fld.scn.env.Ωz], [31, 27])

	TL(r, z) = fld.scn.env.ati.z(r) ≤ z ≤ fld.scn.env.bty.z(r) ? min(100.0, fld.TL(r, z)) : NaN

	DEF_NAME = "TL Grid"
	progress_name(name) = length(name) == 0 ? 
	DEF_NAME : name * ": " * DEF_NAME * " "
	pn = progress_name(fld.scn.name)
	TLgrid = @showprogress 1 pn [TL(r′, z′) for z′ ∈ z, r′ ∈ r]

	heatmap!(
		r, z, TLgrid,
		seriescolor = cgrad(:jet, rev = true),
		colorbar = true,
		colorbar_title = "Transmission Loss (dB)"
	)
	
	plot!(title = "Field: " * fld.scn.name)
end

function plot_oac(oac::OceanAcoustic)
	p = plot_oac()
	plot_oac!(oac)
	plot_oac!()
	return p
end

function plot_oac(oac::OceanAcoustic, args::Any...)
	p = plot_oac()
	plot_oac!(oac, args...)
	plot_oac!()
	return p
end

function save_oac_plot(
	p::Plot,
	file_dir_name::AbstractString...
	)
	savefig(p, plotsdir(file_dir_name...) * ".png")
	nothing
end

function save_oac_plot(
	p::Plot,
	file_dir_name...
	)
	fdn = [String(dir) for dir ∈ file_dir_name]
	save_oac_plot(p, fdn...)
	nothing
end
