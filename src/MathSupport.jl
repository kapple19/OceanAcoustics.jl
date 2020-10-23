export interpolated_function
export bivariate_derivatives
export univariate_interpolation
export bivariate_interpolation

function interpolated_function(x, y)
	itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
	@show x
	@show y
	return itp_fcn(x::Real) = itp(x)
end
function interpolated_function(x, y, z)
	itp = LinearInterpolation((x, y), z, extrapolation_bc = Flat())
	return itp_fcn(x::Real, y::Real) = itp(x, y)
end

function bivariate_derivatives(f::Function)
	∂f_∂x(x, y) = derivative(x -> f(x, y), x)
	∂f_∂y(x, y) = derivative(y -> f(x, y), y)

	∂²f_∂x²(x, y) = derivative(x -> ∂f_∂x(x, y), x)
	∂²f_∂y∂x(x, y) = derivative(y -> ∂f_∂x(x, y), y)
	∂²f_∂x∂y(x, y) = derivative(x -> ∂f_∂y(x, y), x)
	∂²f_∂y²(x, y) = derivative(y -> ∂f_∂y(x, y), y)

	return ∂f_∂x, ∂f_∂y, ∂²f_∂x², ∂²f_∂x∂y, ∂²f_∂y∂x, ∂²f_∂y²
end

# using ForwardDiff: gradient
# function bivariate_derivatives(f::Function)
# 	f′(x) = f(x[1], x[2])
# 	∇f′(x) = gradient(f′, x)
# 	∂f_∂x(x, y) = ∇f′([x, y])[1]
# 	∂f_∂y(x, y) = ∇f′([x, y])[2]

# 	∂f_∂x′(x) = ∂f_∂x(x[1], x[2])
# 	∇∂f_∂x′(x) = gradient(∂f_∂x′, x)
# 	∂²f_∂x²(x, y) = ∇∂f_∂x′([x, y])[1]
# 	∂²f_∂y∂x(x, y) = ∇∂f_∂x′([x, y])[2]

# 	∂f_∂y′(x) = ∂f_∂y(x[1], x[2])
# 	∇∂f_∂y′(x) = gradient(∂f_∂y′, x)
# 	∂²f_∂x∂y(x, y) = ∇∂f_∂y′([x, y])[1]
# 	∂²f_∂y²(x, y) = ∇∂f_∂y′([x, y])[2]

# 	return (
# 		∂f_∂x, ∂f_∂y,
# 		∂²f_∂x²,
# 		∂²f_∂x∂y, ∂²f_∂y∂x,
# 		∂²f_∂y²
# 	)
# end

univariate_interpolation(f::Function) = x::Real -> f(x)
univariate_interpolation(v::Real) = x::Real -> v
function univariate_interpolation(
	x::AbstractVector{Tx},
	y::AbstractVector{Ty}
	) where {Tx <: Real, Ty <: Real}
	return interpolated_function(x, y)
end
univariate_interpolation(t::Tuple) = univariate_interpolation(t...)

"""
Returns a bivariate function.
"""
bivariate_interpolation(f::Function) = (x::Real, y::Real) -> f(x, y)

"""
Returns a bivariate constant function.
"""
bivariate_interpolation(v::Real) = (x::Real, y::Real) -> v
function bivariate_interpolation(
	x::AbstractVector{Tx},
	y::AbstractVector{Ty}
	) where {Tx <: Real, Ty <: Real}
	itp = interpolated_function(x, y)
	return (x::Real, y::Real) -> itp(y)
end

"""
Returns a bivariate function defined by a grid of values. Can be irregularly spaced.
"""
function bivariate_interpolation(
	xVec::AbstractVector{Tx},
	yVec::AbstractVector{Ty},
	zVec::AbstractArray{Tz, 2}
	) where {Tx <: Real, Ty <: Real, Tz <: Real}
	return interpolated_function(xVec, yVec, zVec)
end

function bivariate_interpolation(
	xVec::AbstractVector{R},
	zVecFcns::AbstractVector{F}
	) where {R <: Real, F <: Function}
	function f(x::Real, y::Real)
		z_at_y_vals = [zVecFcns[nFcn](y) for nFcn ∈ eachindex(zVecFcns)]
		itp = interpolated_function(xVec, z_at_y_vals)
		return itp(x)
	end
end

"""
Returns a bivariate function interpolated first in depth, then between the provided ranges, i.e. each `(nx, x) ∈ enumerate(xVec)` has a specified profile output `zTup[nx]::AbstractVector` interpolated linearly with `yTup[nx]::AbstractVector`.

`length(x)Vec`

```julia
x = [0, 1, 2, 5]
y = [
	[0, 1, 2, 4, 6],
	[0, 5, 15, 50],
	[0, 100],
	[0, 10, 20, 40]
]
z = []
	[10, 9, 10, 11, 12],
	[5, 4, 3, 4],
	[20, 25],
	[1, 2, 3, 10]
]

f = bivariate_interpolation(x, y, z)
nx = 2, ny = 3

f(x[nx], y[nx][ny]) == z[nx][ny]
```
"""
function bivariate_interpolation(
	xVec::AbstractVector{R},
	yTup::AbstractVector{V},
	zTup::AbstractVector{V}) where {R <: Real, V <: AbstractVector{R}}

	z_wrt_y_fcns = interpolated_function.(yTup, zTup)
	f = bivariate_interpolation(xVec, z_wrt_y_fcns)
	return f
end
