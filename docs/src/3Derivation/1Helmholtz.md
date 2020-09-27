# Helmholtz Equation
The following theory is taken from Jensen et al [ref].

The `OceanAcoustics.jl` package is an implementation of a numerical solution method of the Helmholtz equation.

```math
\nabla^2 p + \frac{\omega^2}{c(\vec{x})^2} p = -\delta(\vec{x} - \vec{x}_0).
```

We assume the solution can be approximated as a ray series

```math
p(\vec{x}) = e^{i\omega\tau(\vec{x})} \sum_{j = 0}^\infty \frac{A_j(\vec{x})}{(i\omega)^{j}}
```

and substitution then equating orders of magnitude, we obtain the equations

```math
\begin{aligned}
& O(\omega^2) &: \left|\nabla\tau\right|^2 &= \frac{1}{c(\vec{x})} \\
& O(\omega) &: 2\nabla\tau \cdot \nabla A_0 + \left( \nabla^2\tau \right) A_0 &= 0 \\
& O(\omega^{1 - j}) &: 2\nabla\tau \cdot A_j + \left( \nabla^2\tau \right) A_j &= -\nabla^2 A_{j-1} & j &= 1, 2, ...
\end{aligned}
```

The ``O(\omega^2)`` equation for ``\tau(\vec{x})``equation is known as the eikonal equation.

The ``O(\omega^j)``, ``j = 0, 1, 2, ...`` equations for ``A_j`` are known as the transport equations.

The implementations in this package approximate a solution for the eikonal equation and the first transport equation.

## Eikonal Equation
Assuming an azimuthal symmetry in cylindrical coordinates, the eikonal equation can be so expressed as to eliminate the ray time ``\tau(s)`` to obtain

```math
\begin{aligned}
\frac{dr}{ds} &= c\xi(s) & \frac{d\xi}{ds} &= \frac{-1}{c}\frac{dc}{dr}, \\
\frac{dz}{ds} &= c\zeta(s) & \frac{d\zeta}{ds} &= \frac{-1}{c}\frac{dc}{dz}
\end{aligned}
```

for the ray trajectory ``\left(r(s), z(s)\right)``, and the ray tangent vector ``c(s) \left(\xi(s), \zeta(s)\right)``. Note that ``\xi(s)`` and ``\zeta(s)`` are in units of slowness (m/s). Additionally, ``c(r(s), z(s)) = c(s)``.

Our initial conditions are

```math
\begin{aligned}
r(0) &= r_0 & \xi(0) &= \frac{\cos(\theta_0)}{c(0)}, \\
z(0) &= z_0 & \zeta(0) &= \frac{\sin(\theta_0)}{c(0)}.
\end{aligned}
```

Additionally, re-expressing the eikonal equation, this time retaining the ray time ``\tau(s)`` we have

```math
\frac{d\tau}{ds} = \frac{1}{c(s)}.
```

with

```math
\tau(0) = \frac{1}{c(0)}
```

## Transport Equation

## Dynamic Ray Equations

## Gaussian Beams