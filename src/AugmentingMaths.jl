using ForwardDiff: gradient
using Interpolations:
LinearInterpolation,
Flat

function interpolated_function(x, y)
	Itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
	return ItpFcn(x::Real) = Itp(x)
end
function interpolated_function(x, y, z)
	Itp = LinearInterpolation((x, y), z, extrapolation_bc = Flat())
	return ItpFcn(x::Real, y::Real) = Itp(x, y)
end

function bivariate_partial_derivatives(f::Function)
	f_(x) = f(x[1], x[2])
	∇f_(x) = gradient(f_, x)
	∇f(x, y) = ∇f_([x, y])
	∂f_∂x(x, y) = ∇f(x, y)[1]
	∂f_∂y(x, y) = ∇f(x, y)[2]

	∂f_∂x_(x) = ∂f_∂x(x[1], x[2])
	∇∂f_∂x_(x) = gradient(∂f_∂x_, x)
	∇∂f_∂x(x, y) = ∇∂f_∂x_([x, y])

	∂f_∂y_(x) = ∂f_∂y(x[1], x[2])
	∇∂f_∂y_(x) = gradient(∂f_∂x_, x)
	∇∂f_∂y(x, y) = ∇∂f_∂y_([x, y])

	∂²f_∂x²(x, y) = ∇∂f_∂x(x, y)[1]
	∂²f_∂x∂y(x, y) = ∇∂f_∂x(x, y)[2]
	∂²f_∂y²(x, y) = ∇∂f_∂y(x, y)[2]

	return ∂f_∂x, ∂f_∂y, ∂²f_∂x², ∂²f_∂x∂y, ∂²f_∂y²
end
