# Implementation Method
This chapter details methods of implementation for the mathematical theory in the OceanAcoustics package.

```@contents
Pages = ["3implementation.md"]
```

## Boundary

## Medium

## Environment

## Position

## Signal

## Fan

## Source

## Scenario

## Trace

## Beam

## Field

### Ray Coordinates to Cylindrical Coordinates
The solution to the ray tracing provides a coordinates system of the vertical slice of ocean in terms of the ray length $s$ and normal $n$. A conversion to cylindrical coordinates for range $r$ and depth $z$ is necessary to proceed.

Given a point $(r', z')$ in the field, we seek $(s, n)$ for a ray whose coordinates are given by $(r(s), z(s))$.

We define the distance quadrature $Q(s)$ as the square of the distance of the field point from the the ray at ray length $s$,

```math
\newcommand{\Par}[1]{\left({#1}\right)}
Q(s) = \Par{r(s) - r'}^2 + \Par{z(s) - z'}^2
```

and implement automatic differentiation to find its minimum value denoted $s'$. The respective normal is thus $n' = \sqrt{Q(s')}$.
