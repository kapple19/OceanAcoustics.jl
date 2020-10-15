src = OceanAcoustics.Source(OceanAcoustics.Position(r₀, z₀), OceanAcoustics.Signal(f))
ocn = OceanAcoustics.Medium(c, R, Z)
bty = OceanAcoustics.Boundary(Z, R)
ati = OceanAcoustics.Boundary(0.0, R)