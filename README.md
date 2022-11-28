# OceanAcoustics
An implementation of acoustics models in the context of long-range ocean propagation.

![Munk Profile Ray Trace](test/img/trace_munk_profile.png)

## Development Roadmap
* Models:
  * Ray/beam tracing
  * Parabolic equation
  * Sonar equations
* Auxiliary:
  * Example scenarios
  * Plot recipes

### Tasks
* Check initial rays launch within ocean.
* Implement more example `Scenario`s.
  * Balearic sea [Jensen et al. p 171]
* Implement beam tracing, then define plot recipe.
* Implement gridded field approximation, then define plot recipe.
* Document model features, comparing with and referencing literature.
* Outline mathematics - replace TeX doc?