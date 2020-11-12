
parse_interval(Ω::Interval) = Ω
parse_interval(V::Real) = V..V

gridpoints(Ω_lo::Real, Ω_hi::Real, N::Integer = DEFAULT_GRID_N) = LinRange(Ω_lo, Ω_hi, N)

gridpoints(Ω::Tuple{Rlo, Rhi}, N::Integer = DEFAULT_GRID_N) where {Rlo <: Real, Rhi <: Real} = gridpoints(Ω[1], Ω[2], N)

gridpoints(Ω::Interval, N::Integer = DEFAULT_GRID_N) = gridpoints(Ω.lo, Ω.hi, N)

function interpolated_function(x, y)
	Itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
	return ItpFcn(x::Real) = Itp(x)
end
function interpolated_function(x, y, z)
	Itp = LinearInterpolation((x, y), z, extrapolation_bc = Flat())
	return ItpFcn(x::Real, y::Real) = Itp(x, y)
end

function closest_points(x′, y′, x, y, Ωs)
	Q(s) = (x(s) - x′)^2 + (y(s) - y′)^2
	dQ(s) = derivative(Q, s)
	return [
		(s, √Q(s))
		for s ∈ find_zeros(dQ, Ωs.lo, Ωs.hi)
	]
end
