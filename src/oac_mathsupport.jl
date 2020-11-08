
parse_interval(Ω::Interval) = Ω
parse_interval(V::Real) = V..V

function interpolated_function(x, y)
	Itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
	return ItpFcn(x::Real) = Itp(x)
end
function interpolated_function(x, y, z)
	Itp = LinearInterpolation((x, y), z, extrapolation_bc = Flat())
	return ItpFcn(x::Real, y::Real) = Itp(x, y)
end

function closest_points(r, z, beam)
	Q(s) = (beam.ray.r(s) - r)^2 + (beam.ray.z(s) - z)^2
	dQ(s) = derivative(Q, s)
	return [
		(s, √Q(s))
		for s ∈ find_zeros(dQ, beam.ray.Ωs.lo, beam.ray.Ωs.hi)
	]
end
