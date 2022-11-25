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
struct Depth <: OAC
	fcn::Function
	min::Float64
	max::Float64
end

export Depth

function Depth(z::Function, domain::AbstractInterval{<:Real})
	rng = z(domain)
	Depth(z, rng.lo, rng.hi)
end

function Depth(z::Function, domain::Tuple{<:Real, <:Real})
	Depth(z, Interval(domain...))
end

function Depth(fcn::Function, min::Real, max::Real)
	zFcn(x::Real) = fcn(x)
	Depth(zFcn, Float64(min), Float64(max))
end

function Depth(x::Vector{<:Real}, z::Vector{<:Real})
	# Checks
	!issorted(x) && throw(NotSorted(x))
	!allunique(x) && throw(NotAllUnique(x))
	length(x) ≠ length(z) && throw(DimensionMismatch())

	Depth(linear_interp_fcn(x, z), minimum(z), maximum(z))
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
mutable struct Surface <: OAC
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
mutable struct Bottom <: OAC
	z::Depth

	function Bottom(args...)
		z = Depth(args...)
		new(z)
	end
end

export Bottom

# Bottom(btm::Bottom) = btm

"""
`Ocean`
"""
mutable struct Ocean <: OAC
	c::Function
	Ocean(c::Function) = new((x, z) -> c(x, z))
end

function Ocean(c::Real)
	Ocean((x, z) -> c)
end

function Ocean(z::AbstractVector{<:Real}, c::AbstractVector{<:Real})
	Ocean(linear_interp_fcn(z, c))
end

export Ocean

Ocean(ocn::Ocean) = ocn

# Ocean(ocn) = Ocean(ocn...)

"""
`Environment`
"""
mutable struct Environment <: OAC
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
mutable struct Source <: OAC
	f::Float64
	z::Float64
end

export Source

Source(src::Source) = src

Source(src) = Source(src...)

"""
`Receiver`
"""
mutable struct Receiver <: OAC
	x::Float64

	Receiver(x::Float64) = new(x)
end

export Receiver

Receiver(rcv::Receiver) = rcv

"""
`Entities`
"""
mutable struct Entities <: OAC
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
mutable struct Scenario <: OAC
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