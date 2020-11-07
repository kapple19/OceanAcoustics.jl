# Mathematical Background
This section expounds the mathematical background implemented in the OceanAcoustics package.

## Ray Theory
The following theory is taken from Jensen.

We are interested in propagation within a slice of ocean. We take cylindrical coordinates with rotational symmetry, with range $r$ and depth $z$ (positive downwards). The sound speed field is denoted as $c(r, z)$, a function of range and depth. When speaking of the coordinates of a ray, $(r(s), z(s))$ is used, $s$ being the path-length with $s = 0$ at the ray origin. Likewise, univariate sound speed $c(s)$ denotes said parameter along the ray. The angle of the ray at ray length $s$ is denoted by $\theta(s)$. We will also be observing normal distance from the ray $n$.

Acoustic propagation in the ocean is modelled with the **wave equation**. Under a frequency-time Fourier transform, the **Helmholtz equation** is obtained. While there are a number of ocean variables that represent wave motion (e.g. pressure, particle velocity, velocity potential, etc.) we focus on pressure.

Perturbation analysis via expressing the solution as a ray series yields the **eikonal equation**, and an infinite set of **transport equations**. The **dynamic ray equations** are employed to complete the calculation of the ray.

### Eikonal Equation
The **eikonal equation** is re-expressed and parameterised with respect to ray length $s$, and the second-order derivatives of $s$ are split into respective systems of equations in terms of the local ray tangent $c(s) (\xi(s), \zeta(s))$. This yields

```math
\newcommand{\dif}[1]{\textup{d}{#1}}
\newcommand{\Dif}[2]{\dfrac{\dif{#1}}{\dif{#2}}}
\newcommand{\Part}[2]{\dfrac{\partial{#1}}{\partial{#2}}}
\begin{aligned}
	\Dif{r}{s} &= c(s) \xi(s), \\
	\Dif{z}{s} &= c(s) \zeta(s), \\
	\Dif{\xi}{s} &= \frac{-1}{c(s)^2} \Part{c}{r}, \\
	\Dif{\zeta}{s} &= \frac{-1}{c(s)^2} \Part{c}{z}, \\
	\Dif{\tau}{s} &= \frac{1}{c(s)}.
\end{aligned}
```

For a ray with launch angle $\theta_0$ originating from a source position $(r_0, z_0)$ we have initial conditions

```math
\begin{aligned}
	c(0) &= c_0 = c(r_0, z_0), \\
	\theta(0) &= \theta_0, \\
	r(0) &= r_0, \\
	z(0) &= z_0, \\
	\xi(0) &= \xi_0 = \frac{\cos(\theta_0)}{c(0)}, \\
	\zeta(0) &= \zeta_0 = \frac{\sin(\theta_0)}{c(0)}, \\
	\tau(0) &= \tau_0 = 0,
\end{aligned}
```
which now completes an initial-value problem.

### Dynamic Ray Equations
To complete the ray pressure calculation, we introduce the **dynamic ray equations** which enable calculation of the amplitude along the ray.

```math
\newcommand{\dif}[1]{\textup{d}{#1}}
\newcommand{\Dif}[2]{\dfrac{\dif{#1}}{\dif{#2}}}
\begin{aligned}
	\Dif{q}{s} &= c(s) p(s), \\
	\Dif{p}{s} &= \frac{-c_{nn}(s)}{c(s)^2} q(s).
\end{aligned}
```

where $c_{nn}$ is the second derivative of the sound speed in a direction normal to the ray path.

```math
\newcommand{\Par}[1]{\left({#1}\right)}
\newcommand{\Partt}[2]{\dfrac{\partial^2{#1}}{\partial{#2}^2}}
\newcommand{\Partts}[3]{\dfrac{\partial^2{#1}}{\partial{#2}\partial{#3}}}
c_{nn}(s) = c(s)^2 \Par{\zeta^2 \Partt{c}{r} - 2\xi\zeta \Partts{c}{r}{z} + \xi^2 \Partt{c}{z}}.
```

The initial conditions are at the discretion of the modeller.

#### Angle Perturbation
One option for initial conditions is perturbation with respect to the ray angle.

```math
\begin{aligned}
	p(0) &= p_0 = \frac{1}{c_0},
	q(0) &= q_0 = 0
\end{aligned}
```

and the amplitude for the first-order transport equation is thus given by

```math
\newcommand{\Abs}[1]{\left\lvert{#1}\right\rvert}
A(s) = \frac{1}{4\pi} \sqrt{\Abs{\frac{c(s)}{c_0} \frac{\cos(\theta_0)}{r(s) q(s)}}}
```

#### 

## Beam Tracing
The pressure field for a slice of ocean is taken via the summation (with choice of coherence) of Gaussian beams, with each traced ray taken to be the center of a beam. The pressure field contribution of a single beam is thus given by

```math
\newcommand{\Par}[1]{\left({#1}\right)}
\newcommand{\Brace}[1]{\left\{{#1}\right\}}
p^\textup{beam}(s, n) = A(s)\sqrt{\frac{c(s)}{r(s) q(s)}} \exp\Brace{-i\omega\Par{\tau(s) + \frac{p(s)}{q(s)} \frac{n^2}{2}}},
```

with $n = 0$ yielding the pressure along the ray.

The given beam pressure must then be expressed in cylindrical coordinates $(r, z)$.

For coherent summation of $M$ beams,

```math
p^{c} &= \sum_{m = 1}^M p_m(r, z)
```
