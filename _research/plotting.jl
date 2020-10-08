##
using OceanAcoustics
using Plots:
plot,
heatmap,
cgrad,
plot!

##
include("../scripts/scenarios.jl")

θ₀, src, ocn, bty, ati, title = simple()

fld = Field(θ₀, src, ocn, bty, ati)

## Labels
function acoustic_plot!()
	plot!(
		xaxis = "Range (m)",
		yaxis = ("Depth (m)", :flip)
	)
end

function acoustic_plot()
	p = plot()
	acoustic_plot!()
	return p
end

acoustic_plot()

## Limits
function acoustic_plot(R::Real, Z::Real)
	p = acoustic_plot()
	plot!(
		xlims = (0, R),
		ylims = (0, Z)
	)
	return p
end

@show ocn.R
@show ocn.Z
acoustic_plot(ocn.R, ocn.Z)

## Boundary
function acoustic_plot!(rng::AbstractVector{T}, bnd::Boundary) where T <: Real
	plot!(rng, bnd.z,
		linecolor = :black)
end

function acoustic_plot(rng::AbstractVector{T}, bnd::Boundary, Z::Real) where T <: Real
	p = acoustic_plot(rng[end], Z)
	acoustic_plot!(rng, bnd)
	return p
end

rng = range(0, fld.ocn.R, length = 101)
acoustic_plot(rng, ati, ocn.Z)
acoustic_plot!(rng, bty)

## Rays
function acoustic_plot!(ray::Ray)
	plot!(ray.sol, vars = (1, 2))
end

function acoustic_plot!(rays::AbstractVector{T}) where T <: Ray
	acoustic_plot!.(rays)
end

function acoustic_plot(ray::Union{Ray, AbstractVector{T}}) where T <: Ray
	p = acoustic_plot()
	acoustic_plot!(ray)
	return p
end

rays = Vector{Ray}(undef, 0)
for nRay = 1:length(fld.beams)
	push!(rays, fld.beams[nRay].ray)
end
acoustic_plot(rays)

## Field
function acoustic_plot(rng::AbstractVector{T}, dpt::AbstractVector{T}, fld::Field) where T <: Real
	p = heatmap(rng, dpt, fld.TL,
		seriescolor = cgrad(:jet, rev = true),
		legend = false,
		xaxis = ("Range (m)", extrema(rng)),
		yaxis = ("Depth (m)", :flip, extrema(dpt)),
		colorbar = :right)
	acoustic_plot!(rng[end], dpt[end])
	# plot!(rng, fld.ati.z)
	# plot!(rng, fld.bty.z)
	# for nRay = 1:length(beams)
	# 	plot!(fld.beams[nRay].ray.sol, vars = (1, 2))
	# end
	return p
end

rng = range(0, ocn.R, length = 31)
dpt = range(0, ocn.Z, length = 15)

acoustic_plot(rng, dpt, fld)

## Altogether
rng = range(0, ocn.R, length = 31)
dpt = range(0, ocn.Z, length = 15)
rays = Vector{Ray}(undef, 0)
for nRay = 1:length(fld.beams)
	push!(rays, fld.beams[nRay].ray)
end
acoustic_plot(rays)

p = acoustic_plot(rng, dpt, fld)
acoustic_plot!(rng, ati)
acoustic_plot!(rng, bty)
acoustic_plot!(rays)
display(p)