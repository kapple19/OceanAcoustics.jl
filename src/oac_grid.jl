export Propagation

abstract type Grid <: OceanAcoustic end

# include("pressure_sequential.jl")
include("pressure_parallel.jl")

Nr = 31
Nz = 29

struct Propagation <: Grid
	scn::Scenario
	r::AbstractVector{Rr} where Rr <: Real
	z::AbstractVector{Rz} where Rz <: Real
	p::AbstractArray{Cp, 2} where Cp <: Number
	TL::AbstractArray{RTL, 2} where RTL <: Real

	function Propagation(
		fld::Field,
		r::AbstractVector{Rr},
		z::AbstractVector{Rz}
		) where {Rr <: Real, Rz <: Real}
		
		DEF_NAME = "Pressure Propagation"
		progress_name(name) = length(name) == 0 ? DEF_NAME : name * ": " * DEF_NAME * " "
		pn = progress_name(fld.scn.name)
	
		function pressure(r, z)
			if fld.scn.env.ati.z(r) < z < fld.scn.env.bty.z(r)
				return fld.p(r, z)
			else
				return NaN
			end
		end
	
		p = calc_pressure_grid(r, z, pressure)

		TL = min.(90, SonarEqs.transmission_loss.(4π*p))
	
		return new(fld.scn, r, z, p, TL)
	end
end

function Propagation(fld::Field, Nr::Integer, Nz::Integer)
	r = gridpoints(fld.scn.env.Ωr, Nr)
	z = gridpoints(fld.scn.env.Ωz, Nz)

	return Propagation(fld, r, z)
end

function Propagation(scn::Scenario, args...)
	fld = Field(scn)
	return Propagation(fld, args...)
end

Propagation(oac::OceanAcoustic) = Propagation(oac, Nr, Nz)

