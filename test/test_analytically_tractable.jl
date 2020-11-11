using OceanAcoustics

z_c = [0, 3e2, 4e3]
c = [1438, 1460, 1519.2]
Z = z_c[end]
R = 1e6

ocn = Medium(z_c, c)
bty = Boundary(Z)
env = Environment(R, ocn, bty)

r₀, z₀ = 0.0, 0.0
θ₀_crit = ocn.SSP.c(r₀, z₀)/ocn.SSP.c(r₀, Z) |> acos
θ₀s = θ₀_crit

fan = Fan(θ₀s)
src = Source(Position(r₀, z₀), Signal(200), fan)
scn = Scenario(env, src, "Simple Convergence Zone")

trc = Trace(scn)

plot_oac(trc)
