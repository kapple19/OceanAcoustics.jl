using OceanAcoustics

scn = ExampleScenarios.wavy()

# sols = OceanAcoustics.propagation(scn)

trc = Trace(scn)

##
using GRUtils

function plot_oac()
	f = Figure()
	hold!(f, true)
	return f
end

function plot_oac!()
	f = gcf()

	yflip!(f, true)
	xlabel!(f, "Range (m)")
	ylabel!(f, "Depth (m)")
end

function plot_oac!(env::Environment)
	f = gcf()

	r = LinRange(
		env.Ωr.lo,
		env.Ωr.hi,
		1001
	)
	plot!(f, r, env.ati.z.(r))
	plot!(f, r, env.bty.z.(r))
end

function plot_oac(env::Environment)
	f = plot_oac()
	plot_oac!(env)
	plot_oac!()
	return f
end

function plot_oac!(ray::Ray)
	f = gcf()

	s = LinRange(ray.Ωₛ.lo, ray.Ωₛ.hi, 1001)
	plot!(f, ray.r.(s), ray.z.(s))
end

function plot_oac(ray::Ray)
	f = plot_oac()
	plot_oac!(ray)
	plot_oac!()
	return f
end

function plot_oac!(trc::Trace)
	f = gcf()
	plot_oac!(trc.scn.env)
	for nRay ∈ 1:length(trc.rays)
		plot_oac!(trc.rays[nRay])
	end
end

function plot_oac(trc::Trace)
	f = plot_oac()
	plot_oac!(trc)
	plot_oac!()
	return f
end

##
f = plot_oac(trc)

# display(f)

savefig("trace_eg.png", f)
