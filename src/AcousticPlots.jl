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

export oac_plot
export oac_plot!

function oac_plot(
	env::Environment,
	rays::AbstractVector{R}
	) where R <: Ray

	r = LinRange(env.Ωrange.lo, env.Ωrange.hi, 1001)
	z = LinRange(env.Ωdepth.lo, env.Ωdepth.hi, 1001)

	f = Figure()

	contourf!(f, r, z, env.c.(r', z))
	oac_plot!.(rays)
	return f
end

function oac_plot!(ray::Environment)
	f = gcf()
	s = LinRange(0, ray.S, 1001)
	plot!(ray.r(s), ray.z(s))
end

# function oac_plot()
# 	f = Figure()
# 	hold!(f, true)
# 	xlabel!(f, "Range (m)")
# 	ylabel!(f, "Depth (m)")
# 	colorscheme!(f, "light")
# 	return f
# end

# function oac_plot!(title::AbstractString)
# 	f = gcf()
# 	title!(f, title)
# end

# function oac_plot!()
# 	f = gcf()
# 	yflip!(f, true)
# end

# function oac_plot!(bnd::Boundary)
# 	f = gcf()
# 	r = LinRange(0.0, bnd.R, 1001)
# 	plot!(f, r, bnd.z,
# 		linecolor = color(0, 0, 0))
# end

# function oac_plot(bnd::Boundary)
# 	f = oac_plot()
# 	oac_plot!(bnd)
# 	oac_plot!()
# 	return f
# end

# function oac_plot!(env::Environment)
# 	oac_plot!.(env.media)
# 	oac_plot!.(env.boundaries)
# end

# function oac_plot(env::Environment)
# 	f = oac_plot()
# 	oac_plot!(env)
# 	return f
# end

# function oac_plot!(ray::Ray)
# 	f = gcf()
# 	s = LinRange(0.0, ray.S, 1001)
# 	plot!(f, ray.r(s), ray.z(s))
# end

# function oac_plot(ray::Ray)
# 	f = oac_plot()
# 	oac_plot!(ray)
# 	oac_plot!()
# 	return f
# end

# function oac_plot!(fld::Field)
# 	f = gcf()
# 	r = LinRange(0.0, fld.ocn.R, 11)
# 	z = LinRange(0.0, fld.ocn.Z, 5)
# 	contourf(r, z, fld.TL.(r', z),
# 		levels = 21, majorlevels = 2)
# end

# function oac_plot(fld::Field)
# 	f = oac_plot()
# 	oac_plot!(fld)
# 	oac_plot!()
# 	return f
# end
