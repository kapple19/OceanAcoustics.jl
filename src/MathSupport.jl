# module MathSupport
using ForwardDiff: derivative
using Interpolations:
LinearInterpolation,
Flat

export interpolated_function
export bivariate_derivatives
export univariate_interpolation
export bivariate_interpolation

function interpolated_function(x, y)
	itp = LinearInterpolation(x, y, extrapolation_bc = Flat())
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

struct NotSorted <: Exception
	vec::AbstractVector
	msg::Vector{AbstractString}
	function NotSorted(
		vec::AbstractVector)
		msg = ["Vector not sorted: ", "."]
		return new(vec, msg)
	end
end

Base.showerror(io::IO, e::NotSorted) = print(io, e.msg[1], e.vec, e.msg[2])

struct Duplicates <: Exception
	vec::AbstractVector
	msg::Vector{AbstractString}
	function Duplicates(
		vec::AbstractVector)
		msg = ["Vector has duplicates: ", "."]
		return new(vec, msg)
	end
end

Base.showerror(io::IO, e::Duplicates) = print(io, e.msg[1], e.vec, e.msg[2])

function findnearestindices(x::T, xVec::AbstractVector{T}) where T <: Real
	if !issorted(xVec)
		NotSorted(xVec) |> throw
	end
	if !allunique(xVec)
		Duplicates(xVec) |> throw
	end
	if x ∉ xVec[1]..xVec[end]
		ErrorException("Input value not in input vector range.") |> throw
	end
	if x == xVec[end]
		return nx
	end
	for nx ∈ 1:length(xVec) - 1
		if x == xVec[nx]
			return nx
		elseif x ∈ xVec[nx]..xVec[nx+1]
			return nx, nx+1
		end
	end
	ErrorException("Could not find indices.") |> throw
end

univariate_interpolation(f::Function) = x -> f(x)
univariate_interpolation(f::Real) = x -> f
univariate_interpolation(f::AbstractVector{T}, x::AbstractVector{T}) where T <: Real = interpolated_function(x, f)

bivariate_interpolation(f::Function) = (x, y) -> f(x, y)
bivariate_interpolation(f::Real) = (x, y) -> f
function bivariate_interpolation(f::AbstractVector{T}, y::AbstractVector{T}) where T <: Real
	itp = interpolated_function(y, f)
	fFcn(x, y) = itp(y)
	return fFcn
end
bivariate_interpolation(f::AbstractArray{T}, x::AbstractVector{T}, y::AbstractVector{T}) where T <: Real = interpolated_function((x, y), f)
function bivariate_interpolation(
	fTup::Tuple{V},
	xVec::AbstractVector{T},
	yTup::Tuple{V}) where {T <: Real, V <: AbstractVector{T}}

	fFcns_wrt_y = Vector{Function}(undef, 0)
	for nx ∈ eachindex(xVec)
		push!(fFcns_wrt_y, interpolated_function(yTup[nx], fTup[nx]))
	end

	function f(x, y)
		nx = findnearestindices(x, xVec)
		if length(nx) == 1
			return fFcns_wrt_y[nx](y)
		elseif length(nx) == 2
			xVals = [xVec[nx[1]], xVec[nx[2]]]
			fVals = [fFcns_wrt_y[xInd](y) for xInd ∈ nx]
			return f_wrt_x_at_y = interpolated_function(xVals, fVals)
		end
	end
	return f
end
# end