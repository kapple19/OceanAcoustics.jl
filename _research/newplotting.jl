##
using OceanAcoustics
using GRUtils
# using Plots

#
include("../scripts/scenarios.jl")

#
θ₀, src, ocn, bty, ati, title = n2linear()

##
rays = Ray.(θ₀, src, ocn, bty, ati)

##
f = Figure()
hold(true)
for nRay = 1:length(rays)
	s = range(0.0, rays[nRay].S, length = 1001)
	plot!(f, rays[nRay].r(s), rays[nRay].z(s))
end
yflip(true)
# yflip(true)

##
function acoustic_plot!(f::Figure, bnd::Boundary)
	gcf(f)
	r = range(0.0, bnd.R)
	plot!(f, r, bnd.z)
end

function acoustic_plot(bnd::Boundary)
	f = acoustic_plot()
	acoustic_plot!(f, bnd)
end

##
function acoustic_plot()
	f = Figure()
	hold(true)
	xlabel!(f, "Range (m)")
	ylabel!(f, "Depth (m)")
	return f
end

f = acoustic_plot()

function acoustic_plot!(f::Figure, ray::Ray)
	s = range(0.0, ray.S, length = 1001)
	gcf(f)
	plot(ray.r(s), ray.z(s))
end

Base.broadcastable(m::Figure) = Ref(m)

acoustic_plot!.(f, rays)

display(f)

function acoustic_plot(ray::Ray)
	f = acoustic_plot()
	acoustic_plot!(f, ray)
end

f = acoustic_plot(rays[1])

##
f = Figure()
nRay = 2
s = LinRange(0.0, rays[nRay].S, 1001)
plot!(f, rays[nRay].r(s), rays[nRay].z(s))
yflip(true)

##
nRay = 1
s = LinRange(0.0, rays[nRay].S, 1001)
r = rays[nRay].r(s)
z = rays[nRay].z(s)
plot(r, z)
f = gcf()
nRay = 2
s = LinRange(0.0, rays[nRay].S, 1001)
r = rays[nRay].r(s)
z = rays[nRay].z(s)
plot!(f, rays[nRay].r(s), rays[nRay].z(s))

##
# Of course first you have to load the package
using GRUtils
# Example data
x = LinRange(0, 10, 500)
y = sin.(x.^2) .* exp.(-x)
# Making a line plot is as simple as this:
plot(x, y)
# Then hold the plot to add the envelope...
hold(true)
# The envelope is given in two columns,
# plotted as dashed lines ("--") in black color ("k")
plot(exp.(0:10).^-1 .* [1 -1], "--k")
# Now set the Y-axis limits, and annotate the plot
ylim(-0.5, 0.5)
legend("signal", "envelope")
xlabel("X")
ylabel("Y")
title("Example plot")

##
fld = Field(θ₀, src, ocn, bty, ati)

##
rng = range(0, ocn.R, length = 51)
dpt = range(0, ocn.Z, length = 25)

# TL = [fld.TL(r, z) for z ∈ dpt, r ∈ rng]
TL = fld.TL.(rng', dpt)

##
p = GRUtils.heatmap(TL)
yflip(true)
cmap = colormap("jet")

##
r = src.pos.r
z = src.pos.z .- [1e-12, 1e-13, 0.0]
@show fld.TL.(r, z)
@show fld.p.(r, z)
@show abs.(fld.p.(r, z))
@show abs.(fld.p.(r, z))/4π
@show -20log10.(abs.(fld.p.(r, z))/4π)
@show min.(100.0, -20log10.(abs.(fld.p.(r, z))/4π))

##
rays = Ray()