function fcn(TupVecs::NTuple{N, V}) where {N, V <: AbstractVector}
	x = 1
	for nFcn ∈ eachindex(TupVecs)
		y = TupVecs[nFcn][2]
		println(y)
	end
end

function fcn(TupFcns::Tuple{Vararg{F}}) where F <: Function
	x = 1
	for nFcn ∈ eachindex(TupFcns)
		y = TupFcns[nFcn](x)
		println(y)
	end
end

function fcn(VecFcns::AbstractVector{F}) where F<: Function
	x = 1
	for nFcn ∈ eachindex(VecFcns)
		y = VecFcns[nFcn](x)
		println(y)
	end
end

TupVecs = (
	[1, 3, 5],
	[1, 4, 6, 9],
	[1, 5]
)

TupFcns = (
	x -> x^2,
	x -> x^3,
	x -> x^4
)

VecFcns = [
	x -> x^2,
	x -> x^3,
	x -> x^4
]

fcn(TupVecs)
fcn(VecFcns)
fcn(TupFcns)
