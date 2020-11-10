using Base.Threads:
@threads,
SpinLock,
lock,
unlock,
nthreads
using ProgressMeter:
Progress,
update!,
next!

function calc_pressure_grid(
	r::AbstractVector{Rr},
	z::AbstractVector{Rz},
	pressure::Function
	) where {Rr <: Real, Rz <: Real}
	Nr = length(r)
	Nz = length(z)
	Nrz = Nr*Nz
	rGrid = reshape(
		r' .* ones(Nz, Nr),
		(Nrz,)
	)
	zGrid = reshape(
		z .* ones(Nz, Nr),
		(Nrz,)
	)

	desc = "Pressure Grid (parallel, " * string(nthreads()) * " threads): "

	prog = Progress(Nrz, desc = desc)
	update!(prog, 0)

	p′ = zeros(Complex, Nrz)
	l = SpinLock()
	@threads for n ∈ eachindex(p′, rGrid, zGrid)
		@inbounds p′[n] = pressure(rGrid[n], zGrid[n])
		lock(l)
		next!(prog, step = 1)
		unlock(l)
	end
	p = reshape(p′, Nz, Nr)
end