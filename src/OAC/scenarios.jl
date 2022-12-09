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
struct Depth <: Oac
	fcn::Function
	min::Float64
	max::Float64

	function Depth(z_fcn::Function, z_min::Real, z_max::Real)
		!(z_min ≤ z_max) && error("`min` must be less than `max`.")
		new(x -> z_fcn(x), Float64(z_min), Float64(z_max))
	end
end

export Depth

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
mutable struct Surface <: Oac
	z::Depth

	function Surface(args...)
		z = Depth(args...)
		new(z)
	end
end

export Surface

Surface() = Surface(0)

# Surface(srf::Surface) = srf

"""
`Bottom`
"""
mutable struct Bottom <: Oac
	z::Depth

	function Bottom(args...)
		z = Depth(args...)
		new(z)
	end
end

export Bottom

Bottom(btm::Bottom) = btm

"""
`Ocean`
"""
mutable struct Ocean <: Oac
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

export Ocean

Ocean(ocn::Ocean) = ocn

# Ocean(ocn) = Ocean(ocn...)

"""
`Environment`
"""
mutable struct Environment <: Oac
	ocn::Ocean
	btm::Bottom
	srf::Surface

	function Environment(ocn::Ocean, btm::Bottom, srf::Surface = Surface())
		new(ocn, btm, srf)
	end
end

export Environment

function Environment(ocn, btm, srf = 0)
	Environment(Ocean(ocn), Bottom(btm), Surface(srf))
end

Environment(env::Environment) = env

Environment(env) = Environment(env...)

"""
`Source`
"""
mutable struct Source <: Oac
	f::Float64
	z::Float64
end

export Source

Source(src::Source) = src

Source(src) = Source(src...)

"""
`Receiver`
"""
mutable struct Receiver <: Oac
	x::Float64

	Receiver(x::Float64) = new(x)
end

export Receiver

Receiver(rcv::Receiver) = rcv

"""
`Entities`
"""
mutable struct Entities <: Oac
	src::Source
	rcv::Receiver

	function Entities(src, rcv)
		new(Source(src), Receiver(rcv))
	end
end

export Entities

Entities(ent::Entities) = ent

Entities(ent) = Entities(ent...)

"""
`Scenario`
"""
mutable struct Scenario <: Oac
	env::Environment
	ent::Entities
	name::String

	function Scenario(env, ent, name = "")
		new(Environment(env), Entities(ent), name)
	end
end

export Scenario

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

export scenarioplot

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

# export calc_bnd_range

function calc_ocean_depth_range(scn::Scenario)
	return [
		calc_bnd_range(scn, :srf) |> minimum
		calc_bnd_range(scn, :btm) |> maximum
	]
end

mutable struct Field
	r::Vector{Float64}
	z::Vector{Float64}
	TL::Matrix{Float64}
end

export Field