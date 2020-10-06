## Flat Scenario
function flat()
	src = Source(Position(0, 200), Signal(200))
	ocn = Medium(1500, 2000, 500)
	bty = Boundary(500)
	ati = Boundary(0)

	θ₀ = π/4*range(-1, 1, length = 51)
	
	θ₀, src, ocn, bty, ati, "Flat Environment"
end

function smooth()
	# Altimetry
	zAtiMin = -10
	zAtiMax = 50
	zAti(r) = zAtiMin + (zAtiMax - zAtiMin)*(sin(r/1e3) + 1.)/2

	# Bathymetry
	rBtyPeak = 5e3
	zBtyMax = 1e3
	zBtyMin = 8e2
	Aᵣ = (2rBtyPeak/3)^2/log((9zBtyMax - 11zBtyMin)/(10(zBtyMax - zBtyMin)))
	zBty(r) = zBtyMax - (zBtyMax - zBtyMin)*exp(-(r - rBtyPeak)^2/4e5)

	# Ocean
	rOcnMax = 10e3
	cOcnMin = 1500
	cOcnMax = 1600
	cSolve(r) = [1 zAti(r) zAti(r)^2
		1 (zAti(r) + zBty(r))/2 ((zAti(r) + zBty(r))/2)^2
		1 zBty(r) zBty(r)^2]
	cSolved(r) = cSolve(r)\[cOcnMax, cOcnMin, cOcnMax]
	cCoeff₀(r) = cSolved(r)[1]
	cCoeff₁(r) = cSolved(r)[2]
	cCoeff₂(r) = cSolved(r)[3]
	cOcn(r, z) = cCoeff₀(r) + cCoeff₁(r)*z + cCoeff₂(r)*z^2

	# Source
	r₀ = 0.0
	z₀ = (zBty(r₀) + zAti(r₀))/2
	f = 250

	# Rays
	θ₀ = acos(cOcn(r₀, z₀)/cOcnMax).*(-1.5:0.125:1.5)

	src = Source(Position(0, z₀), Signal(f))
	ocn = Medium(cOcn, rOcnMax, zBtyMax)
	bty = Boundary(zBty)
	ati = Boundary(zAti)
	
	θ₀, src, ocn, bty, ati, "Smooth Environment"
end

## Convergence Zones
function convergence()
	c = [1520, 1500, 1515, 1495, 1545.]
	z = [0., 300., 1200., 2e3, 5000.]
	Z = z[end]
	R = 250e3

	src = Source(Position(0, 0), Signal(200))
	ocn = Medium(z, c, R, Z)
	bty = Boundary(5e3)
	ati = Boundary(0)
	
	θ_crit = acos(ocn.c(0, 0)/ocn.c(0, 5e3))
	θ₀ = θ_crit*range(0.01, 1, length = 10)

	return θ₀, src, ocn, bty, ati, "Convergence Zones"
end

## Upward-Refracting Rays
function upward()
	c(r, z) = 1500 + 100z/5e3
	θ₀_crit = acos(c(0, 0)/c(0, 5e3))
	θ₀ = θ₀_crit*range(0.1, 1, length = 10)
	Z = 5e3

	src = Source(Position(0, 0), Signal(200))
	ocn = Medium(c, 1e5, Z)
	bty = Boundary(Z)
	ati = Boundary(0)
	
	return θ₀, src, ocn, bty, ati, "Upward-Refracting Rays"
end

## EOF
nothing