z_c = [0, 3e2, 4e3]
c = [1438, 1460, 1519.2]
Z = z_c[end]
R = 250e3

ocn = Medium(z_c, c)
bty = Boundary(Z)
env = Environment(R, ocn, bty)

r₀, z₀ = 0.0, 0.0
c₀ = ocn.SSP.c(r₀, z₀)
θ₀_crit = c₀/ocn.SSP.c(r₀, Z) |> acos
θ₀s = θ₀_crit * LinRange(0.5, 1, 31)

fan = Fan(θ₀s)
src = Source(Position(r₀, z₀), Signal(200), fan)
scn = Scenario(env, src, "Simple Convergence Zone")

trc = Trace(scn)

plot_oac(trc)

##
θ₁ = c₀/ocn.SSP.c(r₀, z_c[2]) |> acos
dz = diff(z_c)
dc = diff(c)
g = dc./dz
R_CZ = 2c₀ * sin(θ₁) * sum(g)

##
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
# θ₀s = LinRange(θ₀_crit₁, θ₀_crit₂, 11)
θ₀s = θ₀_crit₂ * LinRange(0.8, 1, 11)

fan = Fan(θ₀s)
src = Source(Position(r₀, z₀), Signal(200.), fan)
scn = Scenario(env, src, "Convergence Zone Propagation")

trc = Trace(scn)

plot_oac(trc)

##