using GRUtils:
Figure,
plot!,
contourf!,
hold!,
yflip!,
xlabel!,
ylabel!,
colorscheme!,
color,
title!,
gcf,
display

## Plots
function acoustic_plot()
	f = Figure()
	hold!(f, true)
	xlabel!(f, "Range (m)")
	ylabel!(f, "Depth (m)")
	colorscheme!(f, "light")
	return f
end

function acoustic_plot!(title::AbstractString)
	f = gcf()
	title!(f, title)
end

function acoustic_plot!()
	f = gcf()
	yflip!(f, true)
end

function acoustic_plot!(bnd::Boundary)
	f = gcf()
	r = LinRange(0.0, bnd.R, 1001)
	plot!(f, r, bnd.z,
		linecolor = color(0, 0, 0))
end

function acoustic_plot(bnd::Boundary)
	f = acoustic_plot()
	acoustic_plot!(bnd)
	acoustic_plot!()
	return f
end

function acoustic_plot!(ray::Ray)
	f = gcf()
	s = LinRange(0.0, ray.S, 1001)
	plot!(f, ray.r(s), ray.z(s))
end

function acoustic_plot(ray::Ray)
	f = acoustic_plot()
	acoustic_plot!(ray)
	acoustic_plot!()
	return f
end

function acoustic_plot!(fld::Field)
	f = gcf()
	r = LinRange(0.0, fld.ocn.R, 11)
	z = LinRange(0.0, fld.ocn.Z, 5)
	contourf(r, z, fld.TL.(r', z),
		levels = 21, majorlevels = 2)
end

function acoustic_plot(fld::Field)
	f = acoustic_plot()
	acoustic_plot!(fld)
	acoustic_plot!()
	return f
end
