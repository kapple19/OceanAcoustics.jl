using ProgressMeter:
@showprogress

function calc_pressure_grid(
	r::AbstractVector{Rr},
	z::AbstractVector{Rz},
	pressure::Function
	) where {Rr <: Real, Rz <: Real}
	p = @showprogress 1 "Pressure Grid (Sequential): " [pressure(r′, z′) for z′ ∈ z, r′ ∈ r]
end
