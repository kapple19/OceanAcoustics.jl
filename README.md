# OceanAcoustics
An implementation of acoustics models in the context of long-range ocean propagation.

![Munk Profile Ray Trace](test/img/trace_munk_profile.png)

## Development Roadmap (Epics)
* Models:
  * Ray/beam tracing
  * Parabolic equation
  * Sonar equations
* Auxiliary:
  * Example scenarios
  * Plot recipes

### Tasks
Dependency-ordered tasks:
* Check initial rays launch within ocean.
* Implement beam tracing, then define plot recipe.
* Implement gridded field approximation, then define plot recipe.
  * Coherent and Incoherent beam summations.
  * Gridded approximation.
* Document model features, comparing with and referencing literature.
  * Tasks will start tree-diverging from here, depending on what features I'm implementing. E.g. boundary losses.
* Outline mathematics - replace TeX doc?
  * LaTeX document: Mathematics.
  * Julia document: Implementation.

Side tasks:
* Docstrings and comments in the source code.
* Colour rays by launch angle magnitude and sign.
* Ribbon the bathymetry and surface plotting.
* Implement more example `Scenario`s.
  * Balearic sea [Jensen et al. p 171]
* Document sonar equations.