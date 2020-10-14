## Flat Scenario
function flat()
	R = 5e2

	src = Source(Position(0.0, 2e2), Signal(2e2))
	ocn = Medium(15e2, 2e3, R)
	bty = Boundary(5e2, R)
	ati = Boundary(0.0, R)

	θ₀ = π/4*LinRange(-1.0, 1.0, 51)
	
	θ₀, src, ocn, bty, ati, "Flat Environment"
end

"""
TODO: Change name to "wavy".
"""
function smooth()
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
	cSolve(r) = [1.0 zAti(r) zAti(r)^2
		1.0 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2.0)^2
		1.0 zBty(r) zBty(r)^2]
	cSolved(r) = cSolve(r)\[cOcnMax, cOcnMin, cOcnMax]
	cCoeff₀(r) = cSolved(r)[1]
	cCoeff₁(r) = cSolved(r)[2]
	cCoeff₂(r) = cSolved(r)[3]
	cOcn(r, z) = cCoeff₀(r) + cCoeff₁(r)*z + cCoeff₂(r)*z^2

	# Source
	r₀ = 0.0
	z₀ = (zBty(r₀) + zAti(r₀))/2
	f = 250.

	# Rays
	θ₀ = acos(cOcn(r₀, z₀)/cOcnMax).*(-1.5:0.125:1.5)

	src = Source(Position(0, z₀), Signal(f))
	ocn = Medium(cOcn, rOcnMax, zBtyMax)
	bty = Boundary(zBty, rOcnMax)
	ati = Boundary(zAti, rOcnMax)
	
	θ₀, src, ocn, bty, ati, "Smooth Environment"
end

## Convergence Zones
function convergence()
	c = [1520, 1500, 1515, 1495, 1545.]
	z = [0., 300., 1200., 2e3, 5000.]
	Z = z[end]
	R = 250e3

	src = Source(Position(0., 0.), Signal(200.))
	ocn = Medium(z, c, R, Z)
	bty = Boundary(5e3, R)
	ati = Boundary(0., R)
	
	θ_crit = acos(ocn.c(0.0, 0.0)/ocn.c(0.0, 5e3))
	θ₀ = θ_crit*LinRange(0.5, 1.0, 10)

	return θ₀, src, ocn, bty, ati, "Convergence Zones"
end

## Upward-Refracting Rays
function upward()
	c(r, z) = 1500 + 100z/5e3
	θ₀_crit = acos(c(0, 0)/c(0, 5e3))
	θ₀ = θ₀_crit*LinRange(0.1, 1.0, 10)
	R = 1e5
	Z = 5e3

	src = Source(Position(0, 0), Signal(200))
	ocn = Medium(c, R, Z)
	bty = Boundary(Z, R)
	ati = Boundary(0, R)
	
	return θ₀, src, ocn, bty, ati, "Upward-Refracting Rays"
end

## Parabolic Bathymetry
function parabolic()
	c = 250.0
	zBty(r) = 2e-3*2.5e5sqrt(1 + r/c)
	θ₀ = LinRange(atan(5e3/2e3), atan(5e3/20e3), 30)
	R = 20e3
	Z = 5e3

	src = Source(Position(0, 0), Signal(200))
	ocn = Medium(c, R, Z)
	bty = Boundary(zBty, R)
	ati = Boundary(0.0, R)

	return θ₀, src, ocn, bty, ati, "Parabolic Bathymetry"
end

## Deep-Sound-Channel
function channel()
	z = [0.0, 500/3, 500/2, 500, 1000, 1500, 4e3]
	c = [1480, 1500, 1485, 1475, 1480, 1485, 1525]
	R = 250e3
	Z = z[end]
	
	z₀ = 500.
	src = Source(Position(0.0, z₀), Signal(200))
	ocn = Medium(z, c, R, Z)
	bty = Boundary(4e3, R)
	ati = Boundary(0.0, R)

	θ₀ = acos(ocn.c(0.0, z₀)/1500.0) * LinRange(-1, 1, 31)

	return θ₀, src, ocn, bty, ati, "Deep Sound Channel"
end

## Seamount
function seamount()
	z = [0, 100, 200, 350, 500, 1500, 3100.]
	c = [1480, 1470, 1475, 1473, 1475, 1488, 1505]
	Z = z[end]
	r = 1e3*[0, 40, 45, 50, 55, 60, 70, 140]
	zBty = [Z, Z, 2900, 2850, 2000, 500, Z, Z]
	R = r[end]

	src = Source(Position(0, 363), Signal(200))
	ocn = Medium(z, c, R, Z)
	bty = Boundary(r, zBty)
	ati = Boundary(0, R)

	θ₀ = atan(363/2e3) * LinRange(-1, 1, 31)

	return θ₀, src, ocn, bty, ati, "Seamount"
end

## Simple Test Scenario
function simple()
	R = 2e3
	Z = 500.
	z = [0, 0.5, 0.6, 0.7, 1]*Z
	c = [1510, 1480, 1500, 1490, 1520]
	rBty = [0, 0.6, 0.7, 0.8, 1]*R
	zBty = [1, 1, 0.8, 1, 1]*Z
	rAti = [0, 0.4, 0.5, 0.6, 1]*R
	zAti = [0, 0, 0.25, 0, 0]*Z
	r₀ = 0.
	z₀ = Z/20.

	src = Source(Position(r₀, z₀), Signal(50.))
	ocn = Medium(z, c, R, Z)
	bty = Boundary(rBty, zBty)
	ati = Boundary(rAti, zAti)

	θ₀_crit = acos(ocn.c(r₀, z₀)/ocn.c(r₀, Z))
	θ₀ = θ₀_crit * (-1.5:0.5:3)

	return θ₀, src, ocn, bty, ati, "Simple Environment"
end

## n²-Linear Profile
function n2linear()
	R = 10e3
	Z = 1e3
	r₀ = 0.0
	z₀ = Z
	f = 2e3
	c₀ = 1550.
	c(r, z) = c₀/sqrt(1 + 2.4z/c₀)

	src = Source(Position(r₀, z₀), Signal(f))
	ocn = Medium(c, R, Z)
	bty = Boundary(Z, R)
	ati = Boundary(0.0, R)

	θ₀_crit = -acos(ocn.c(r₀, z₀)/ocn.c(r₀, 150.))
	θ₀ = θ₀_crit*(0.1:0.05:1.1)

	return θ₀, src, ocn, bty, ati, "n^2-Linear Profile"
end

## Sloped Profiles
function slopes()
	R = 10e3
	Z = 2e3
	r₀ = 0
	z₀ = Z/4
	c(r, z) = 1500 - 100r/R + 100z/Z
	zBty(r) = Z - 500r/R
	zAti(r) = 100r/R

	src = Source(Position(r₀, z₀), Signal(200))
	ocn = Medium(c, R, Z)
	bty = Boundary(zBty, R)
	ati = Boundary(zAti, R)

	θ₀ = acos(c(r₀, z₀)/c(r₀, Z)) * (-2:0.2:2)
	
	return θ₀, src, ocn, bty, ati, "Linearly Sloped Profiles"
end

## EOF
nothing