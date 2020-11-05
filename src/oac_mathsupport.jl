
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
	sMins = find_zeros(dQ, beam.ray.Ωs.lo, beam.ray.Ωs.hi)
	d²Q(s) = derivative(dQ, s)
	# min_cond(s) = d²Q(s) > 0 && beam.W(s) > sqrt(Q(s))
	min_cond(s) = d²Q(s) > 0
	min_cond.(sMins)
	filter!(min_cond, sMins)
	# return sMins, sqrt.(Q.(sMins))
	return [(s, √Q(s)) for s ∈ sMins]
end
