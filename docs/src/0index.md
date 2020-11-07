# Home
This documentation details the implementations of ocean acoustics modelling in the package OceanAcoustics.jl.

## Documentation Outline
```@contents
Depth = 1
```

## Summary of Implementation Features
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

Performance tools:
* Multithreading, parallelism on CPU/GPU

User-friendly experience:
* Simple and flexible code use
* Progress bar for long calculations (ray tracing, transmission loss)
