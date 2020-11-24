##
using OceanAcoustics
name = :channel
scn = example_scenario(name)
trc = Trace(scn)
ray = trc.rays[1];
beam = Beam(ray, scn.src)

##
using NLsolve
using ForwardDiff: derivative

r(s) = ray.r(s)
z(s) = ray.z(s)
dr(s) = derivative(r, s)
dz(s) = derivative(z, s)

function f!(F, x, x₀)
	s, n = x
	r₀, z₀ = x₀
	F[1] = r(s) + n*dz(s) - r₀
	F[2] = z(s) - n*dr(s) - z₀
end

@time nlsolve(f!, [])