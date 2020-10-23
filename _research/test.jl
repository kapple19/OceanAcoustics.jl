## Seamount
using OceanAcoustics
using Plots

##
zc = [0, 100, 200, 350, 500, 1500, 3100.]
c = [1480, 1470, 1475, 1473, 1475, 1488, 1505.]

plot(c, zc, yaxis = :flip)

##
Z = zc[end]

##
rBty = 1e3*[0, 40, 45, 50, 55, 60, 70, 140]
zBty = [Z, Z, 2900, 2850, 2000, 500, Z, Z]

R = rBty[end]

plot(rBty, zBty, yaxis = :flip)

##
ocn = Medium((zc, c))

r = LinRange(0, R, 101)
z = LinRange(0, Z, 101)

heatmap(r, z, ocn.SSP.c, yaxis = :flip)

##
bty = Boundary((rBty, zBty), 1600)

plot(r, bty.z, yaxis = :flip)

##
env = Environment(R, ocn, bty)

Ωr = 0..140000
Ωr |> bty.z

## Scenario
fan = Fan(atan(363/2e3) * LinRange(-1, 1, 31))
src = Source(Position(0, 363), Signal(200), fan)
scn = Scenario(env, src, "Seamount")