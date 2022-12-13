"""
`Depth`

Storing univariate real-valued function (F64 -> F64) and meta-information.

`Depth(z::Function)`
`Depth(z::Function, min::Real, max::Real)`
`Depth(x::Vector{<:Real}, z::Vector{<:Real})` creates an interpolator
`Depth(z::Real)` creates a function

Used for:
	* ocean surface altimetry
	* ocean bottom bathymetry

Author Note: May deprecate. The min/max values storage is used for plotting, but these calculations can be done then instead.
"""
struct Depth <: OACBase.Oac
	fcn::Function
	min::Float64
	max::Float64

	function Depth(z_fcn::Function, z_min::Real, z_max::Real)
		!(z_min ≤ z_max) && error("`min` must be less than `max`.")
		new(x -> z_fcn(x), Float64(z_min), Float64(z_max))
	end
end

function Depth(z::Function, domain::AbstractInterval{<:Real})
	rng = z(domain)
	Depth(z, rng.lo, rng.hi)
end

function Depth(z::Function, domain::Tuple{<:Real, <:Real})
	Depth(z, Interval(domain...))
end

function Depth(x::Vector{<:Real}, z::Vector{<:Real})
	# Checks
	!issorted(x) && throw(NotSorted(x))
	!allunique(x) && throw(NotAllUnique(x))
	length(x) ≠ length(z) && throw(DimensionMismatch())

	z_interp = linear_interpolation(x, z, extrapolation_bc = Line())
	z_fcn(x) = z_interp(x)

	Depth(z_fcn, minimum(z), maximum(z))
end

function Depth(z::Real)
	zF64 = Float64(z)
	Depth(x -> zF64, zF64, zF64)
end

(z::Depth)(x) = z.fcn(x)

# Depth(dpt::Depth) = dpt # not sure this one is needed

# Depth(dpt) = Depth(dpt...) # not sure this one is needed

"""
`Surface`
"""
mutable struct Surface <: OACBase.Oac
	z::Depth

	function Surface(args...)
		z = Depth(args...)
		new(z)
	end
end

Surface() = Surface(0)

# Surface(srf::Surface) = srf

"""
`Bottom`
"""
mutable struct Bottom <: OACBase.Oac
	z::Depth

	function Bottom(args...)
		z = Depth(args...)
		new(z)
	end
end

Bottom(btm::Bottom) = btm

"""
`Ocean`
"""
mutable struct Ocean <: OACBase.Oac
	c::Function
	Ocean(c::Function) = new((x, z) -> c(x, z))
end

function Ocean(c::Real)
	Ocean((x, z) -> c)
end

function Ocean(z::AbstractVector{<:Real}, c::AbstractVector{<:Real})
	c_interp = linear_interpolation(z, c, extrapolation_bc = Line())
	c_fcn(x, z) = c_interp(z)
	Ocean(c_fcn)
end

Ocean(ocn::Ocean) = ocn

# Ocean(ocn) = Ocean(ocn...)

"""
`Environment`
"""
mutable struct Environment <: OACBase.Oac
	ocn::Ocean
	btm::Bottom
	srf::Surface

	function Environment(ocn::Ocean, btm::Bottom, srf::Surface = Surface())
		new(ocn, btm, srf)
	end
end

function Environment(ocn, btm, srf = 0)
	Environment(Ocean(ocn), Bottom(btm), Surface(srf))
end

Environment(env::Environment) = env

Environment(env) = Environment(env...)

"""
`Source`
"""
mutable struct Source <: OACBase.Oac
	f::Float64
	z::Float64
end

Source(src::Source) = src

Source(src) = Source(src...)

"""
`Receiver`
"""
mutable struct Receiver <: OACBase.Oac
	x::Float64

	Receiver(x::Float64) = new(x)
end

Receiver(rcv::Receiver) = rcv

"""
`Entities`
"""
mutable struct Entities <: OACBase.Oac
	src::Source
	rcv::Receiver

	function Entities(src, rcv)
		new(Source(src), Receiver(rcv))
	end
end

Entities(ent::Entities) = ent

Entities(ent) = Entities(ent...)

"""
`Scenario`
"""
mutable struct Scenario <: OACBase.Oac
	env::Environment
	ent::Entities
	name::String

	function Scenario(env, ent, name = "")
		new(Environment(env), Entities(ent), name)
	end
end

Scenario(scn::Scenario) = scn

Scenario(scn) = Scenario(scn...)

@userplot ScenarioPlot
@recipe function plot(sp::ScenarioPlot)
	# Inputs
	scn = sp.args[1]

	# Plot Design
	legend --> :none
	yflip := true
	
	# Extrema
	ocn_depth_range = calc_ocean_depth_range(scn)
	ex = (
		srf = minimum(ocn_depth_range),
		btm = maximum(ocn_depth_range)
	)
	ylims --> (ex.srf, ex.btm)

	# Plot Labels
	plot_title := scn.name
	xguide := "Range [m]"
	yguide := "Depth [m]"

	x = range(0.0, scn.ent.rcv.x)
	# Boundaries
	for boundary in (:srf, :btm)
		bnd = getproperty(scn.env, boundary)
		x = range(0.0, scn.ent.rcv.x)
		z = bnd.z.(x)
		@series begin
			linecolor := :brown
			fillrange := zeros(size(z)) .+ ex[boundary]
			fillcolor := :brown
			x, z
		end
	end
end

function calc_bnd_range(scn::Scenario, bnd::Symbol)
	x_rng = 0.0 .. scn.ent.rcv.x
	z_rng = getproperty(scn.env, bnd).z(x_rng)
	z_rng_int = if !(z_rng isa Interval)
		Interval(z_rng, z_rng)
	else
		z_rng
	end
	z_rng_int.lo, z_rng_int.hi
end

function calc_ocean_depth_range(scn::Scenario)
	return [
		calc_bnd_range(scn, :srf) |> minimum
		calc_bnd_range(scn, :btm) |> maximum
	]
end

mutable struct Field
	r::Vector{Float64}
	z::Vector{Float64}
	PL::Matrix{Float64}
end

@userplot PropagationPlot
@recipe function plot(pp::PropagationPlot)
	fld = pp.args[1]

	legend --> :none
	yflip := true

	@series begin
		seriestype := :heatmap
		seriescolor := cgrad(:jet, rev = true)
		fld.r, fld.z, fld.PL'
	end

end

# Exports
## Types
export Depth
export Surface
export Bottom
export Ocean
export Environment
export Source
export Receiver
export Entities
export Scenario
export Field

## Auxiliary Functions
export calc_ocean_depth_range

## Plot Recipes
export scenarioplot
export scenarioplot!
export propagationplot