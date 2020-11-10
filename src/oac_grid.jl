
# include("pressure_sequential.jl")
include("pressure_parallel.jl")

struct Grid <: OceanAcoustic
	scn::Scenario
	r::AbstractVector{Rr} where Rr <: Real
	z::AbstractVector{Rz} where Rz <: Real
	p::AbstractArray{Cp, 2} where Cp <: Number
	TL::AbstractArray{RTL, 2} where RTL <: Real

	function Grid(
		fld::Field,
		r::AbstractVector{Rr},
		z::AbstractVector{Rz}
		) where {Rr <: Real, Rz <: Real}
		
		DEF_NAME = "Pressure Grid"
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

function Grid(fld::Field, Nr::Integer, Nz::Integer)
	r = gridpoints(fld.scn.env.Ωr, Nr)
	z = gridpoints(fld.scn.env.Ωz, Nz)

	return Grid(fld, r, z)
end

function Grid(scn::Scenario, args...)
	fld = Field(scn)
	return Grid(fld, args...)
end

Grid(oac::OceanAcoustic) = Grid(oac, 31, 27)
