module ExampleScenarios
using OceanAcoustics

function convergence(θ₀_mult = LinRange(0.5, 1.0, 10))
	c = [1522, 1501, 1514, 1496, 1545.]
	z = [0., 300., 1200., 2e3, 5000.]
	Z = z[end]
	R = 250e3

	ocn = Medium(z, c)
	bty = Boundary(5e3)
	ati = Boundary(0.)
	env = Environment(R, ocn, bty, ati)

	r₀ = 0.0
	z₀ = 20.0
	θ₀_crit₂ = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, Z) |> acos
	θ₀_crit₁ = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, 4e3) |> acos
	θ₀s = LinRange(θ₀_crit₁, θ₀_crit₂, 11)

	fan = Fan(θ₀s)
	src = Source(Position(0., 0.), Signal(200.), fan)
	scn = Scenario(env, src, "Convergence Zone Propagation")
end

function wavy(θ₀_mult = -1.5:0.125:1.5)
	# Altimetry
	zAtiMin = 0.
	zAtiMax = 50.
	zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2

	# Bathymetry
	rBtyPeak = 5e3
	zBtyMax = 1e3
	zBtyMin = 8e2
	Aᵣ = (2rBtyPeak/3)^2/log((9.0zBtyMax - 11.0zBtyMin)/(10.0(zBtyMax - zBtyMin)))
	zBty(r) = zBtyMax - (zBtyMax - zBtyMin)*exp(-(r - rBtyPeak)^2/4e5)

	# Ocean
	rOcnMax = 10e3
	cOcnMin = 1500.
	cOcnMax = 1600.
	cSolve(r) = [
		1.0 zAti(r) zAti(r)^2
		1.0 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2)^2
		1.0 zBty(r) zBty(r)^2
	]
	cCoeff(r) = cSolve(r)\[cOcnMax, cOcnMin, cOcnMax]
	cOcn(r, z) = cCoeff(r)' * [1, z, z^2]

	# Environment
	ocn = Medium(cOcn)
	bty = Boundary(zBty)
	ati = Boundary(zAti)
	env = Environment(rOcnMax, ocn, bty, ati)

	# Scenario
	r₀ = 0.0
	z₀ = (zBty(r₀) + zAti(r₀))/2
	θ₀_crit = acos(cOcn(r₀, z₀)/cOcnMax)

	fan = Fan(θ₀_crit * θ₀_mult)
	src = Source(Position(r₀, z₀), Signal(250.), fan)
	scn = Scenario(env, src, "Wavy Environment")
end

function flat(θ₀_mult = LinRange(-1, 1, 51))
	# Environment
	ocn = Medium(1500)
	bty = Boundary(5e2)
	env = Environment(10e2, ocn, bty)

	# Scenario
	fan = Fan(π/4 * θ₀_mult)
	src = Source(Position(0, 2e2), Signal(50), fan)
	scn = Scenario(env, src, "Flat Environment")
end

function slopes(θ₀_mult = -2:0.2:2)
	# Environment
	R = 10e3
	Z = 2e3
	c(r, z) = 1500 - 100r/R + 100z/Z
	zBty(r) = Z - 500r/R
	zAti(r) = 100r/R

	ocn = Medium(c)
	bty = Boundary(zBty)
	ati = Boundary(zAti)
	env = Environment(R, ocn, bty, ati)

	# Scenario
	r₀ = 0
	z₀ = Z/4

	fan = Fan(acos(c(r₀, z₀)/c(r₀, Z)) * θ₀_mult)
	src = Source(Position(r₀, z₀), Signal(50), fan)
	scn = Scenario(env, src, "Sloped Environment")
end

function parabolic()
	# Environment
	c = 250
	b = 2.5e5
	zBty(r) = 2e-3b * √(1 + r/c)
	R = 20e3
	Z = 5e3

	ocn = Medium(c)
	bty = Boundary(zBty)
	env = Environment(R, ocn, bty)

	# Scenario
	fan = Fan(LinRange(atan(5e3/2e3), atan(5e3/26e3), 31))
	src = Source(Position(0, 0), Signal(50), fan)
	scn = Scenario(env, src, "Parabolic Bathymetry")
end

function upward(θ₀_mult = LinRange(0.1, 1.0, 10))
	# Environment
	R = 1e5
	Z = 5e3
	c(r, z) = 1500 + 100z/Z

	ocn = Medium(c)
	bty = Boundary(Z)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 0

	fan = Fan(acos(c(r₀, z₀)/c(r₀, Z)) * θ₀_mult)
	src = Source(Position(r₀, z₀), Signal(100), fan)
	scn = Scenario(env, src, "Upward-Refracting Rays")
end

function seamount(θ₀_mult = LinRange(-1, 1, 31))
	# Environment
	zc = [0, 100, 200, 350, 500, 1500, 3100.]
	c = [1480, 1470, 1475, 1473, 1475, 1488, 1505.]

	Z = zc[end]

	rBty = 1e3*[0, 40, 45, 50, 55, 60, 70, 140]
	zBty = [Z, Z, 2900, 2850, 2000, 500, Z, Z]

	R = rBty[end]

	ocn = Medium(zc, c)
	bty = Boundary(rBty, zBty)
	env = Environment(R, ocn, bty)

	# Scenario
	fan = Fan(atan(363/2e3) * θ₀_mult)
	src = Source(Position(0, 363), Signal(200), fan)
	scn = Scenario(env, src, "Seamount")
end

function channel(θ₀_mult = LinRange(-1, 1, 31))
	# Environment
	z = [0.0, 500/3, 500/2, 500, 1000, 1500, 4e3]
	c = [1480, 1500, 1485, 1475, 1480, 1485, 1525.]
	R = 250e3
	Z = z[end]

	ocn = Medium(z, c)
	bty = Boundary(Z)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 500
	θ₀_crit = ocn.SSP.c(r₀, z₀)/1500.0 |> acos

	fan = Fan(θ₀_crit * θ₀_mult)
	src = Source(Position(r₀, z₀), Signal(150), fan)
	scn = Scenario(env, src, "Deep Sound Channel")
end

function munk(θ₀_mult = LinRange(-1, 1, 31))
	# Environment
	z_(z) = 2(z - 1300)/1300
	ϵ = 7.37e-3
	c(r, z) = 1500*(1 + ϵ*(z_(z) - 1 + exp(-z_(z))))

	ocn = Medium(c)
	bty = Boundary(5e3)
	env = Environment(100e3, ocn, bty)

	# Scenario
	r₀ = 0.0
	z₀ = 1e3
	θ₀_crit = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, 0) |> acos

	fan = Fan(θ₀_crit * θ₀_mult)
	src = Source(Position(r₀, z₀), Signal(150), fan)
	scn = Scenario(env, src, "Munk Profile")
end

function n2linear(θ₀_mult = LinRange(0.8, 1.2, 21))
	# Environment
	c₀ = 1550
	c(r, z) = c₀/√(1 + 2.4z/c₀)
	R = 4e3
	Z = 1e3

	ocn = Medium(c)
	bty = Boundary(Z)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0.0
	z₀ = Z
	θ₀_crit = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, 0) |> acos
	θ₀s = -θ₀_crit * θ₀_mult

	fan = Fan(θ₀s)
	src = Source(Position(r₀, z₀), Signal(2e3), fan)
	scn = Scenario(env, src, "n²-Linear Profile")
end

function stacked(θ₀_mult = LinRange(0, 1, 11))
	@warn "Stacked layers not yet supported."

	# Environment
	Z = 5e3
	zStacks = Z*(0:0.2:1)
	cStacks = 1470:10:1530
	function c(r, z)
		for (nz, z′) ∈ enumerate(zStacks)
			if z ≤ z′
				return cStacks[nz]
			end
		end
		return cStacks[end]
	end
	R = 10e3

	ocn = Medium(c)
	bty = Boundary(Z)
	env = Environment(R, ocn, bty)

	# Scenario
	r₀ = 0
	z₀ = 0
	θ_bot = π/4
	θ₀_crit = c(r₀, z₀)/c(r₀, Z)*cos(θ_bot) |> acos
	θ₀s = θ₀_crit * θ₀_mult

	fan = Fan(θ₀s)
	src = Source(Position(r₀, z₀), Signal(50), fan)
	scn = Scenario(env, src, "Stacked Layers")
end

# Export all (at end of module)
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end