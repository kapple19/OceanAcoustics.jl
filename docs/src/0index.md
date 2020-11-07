# Home
This documentation details the implementations of ocean acoustics modelling in the package OceanAcoustics.jl.

## Documentation Outline
```@contents
Depth = 1
```

## Installation
To install this package, simple enter
```julia
add OceanAcoustics
```
in Julia's `Pkg` mode, accessed by pressing `;` from Julian mode.

## Load
To load this package, simple enter
```julia
using OceanAcoustics
```
in any file or in Julian mode.

## Simple Demonstration of Features

```@example
using OceanAcoustics

Z = 5e3
R = 10e3
bty_z(r) = Z - Z/2.5*exp(-(r - R/2)^2/1e5)
ati_z(r) = 200r/R
c(r, z) = 1500 + 100z/Z - 200r/R

bty = Boundary(bty_z)
ati = Boundary(ati_z)
ocn = Medium(c)
env = Environment(R, ocn, bty, ati)

r₀ = 0
z₀ = 3e3
θ_crit = c(r₀, z₀)/c(r₀, Z) |> acos

pos = Position(r₀, z₀)
fan = LinRange(-π/3, π/3, 21) |> Fan
sig = Signal(100)
src = Source(pos, sig, fan)
scn = Scenario(env, src, "Documentation Example")

trc = Trace(scn)

plot_oac(trc)
```

An acoustic ray method in a 2D slice of ocean:
* Range-dependent ocean boundaries
* Range- and depth-dependent sound speed profile
* Ray tracing

Acoustic slice field computations:
* Beam tracing
* Pressure field
* Transmission loss

Detection theory:
* Detection index
* Detection threshold
* Signal excess
* Probability of detection

Built-in examples:
* Scenarios
* Traces
* Beams
* Sonar

Performance tools:
* Multithreading, parallelism on CPU/GPU

User-friendly experience:
* Simple and flexible code use
* Progress bar for long calculations (ray tracing, transmission loss)
