# Ray Tracing
The environment modelled is a vertical slice of ocean with range and depth dependent sound speed profile, bathymetry, and altimetry. Acoustic rays propagate from a single point source. All dynamics are static.

## Environment
The ocean environment is defined with Julia `struct`s.

### Medium
The volume of ocean within which rays propagate is termed a `Medium`. This struct receives two medium properties:
* Range
* Sound speed profile

The range is implemented to abort the ray tracing upon reaching this specified range.

The sound speed profile can be defined as either:
* A single value
* A range-independent, depth-dependent pair of vectors for depth and sound speed profile
* Range and depth vectors accompanying a matrix sound speed profile
* A bivariate function of range and depth

As a function:

```@example
using OceanAcoustics # hide
R = 5e3 # range (m)
c(r, z) = 1500 + 1e-2z^2 # sound speed (m/s)
ocn = Medium(c, R)
```

As a constant:

```@example
using OceanAcoustics # hide
R = 10e3
c = 1520
ocn = Medium(c, R)
```

### Boundary
The volume of ocean is bounded above by the altimetry, and below by the bathymetry. Each are represented as a `Boundary`.

```@example
using OceanAcoustics # hide
ati = Boundary(0)
bty = Boundary(1e3)
```

### Signal
The acoustic source is partially defined by the properties of the source `Signal`.

```@example
using OceanAcoustics # hide
f = 1e3 # frequency (Hz)
sig = Signal(f)
```

### Position
The position of an object in the ocean is 

### Sound Source


## Ray Trace

