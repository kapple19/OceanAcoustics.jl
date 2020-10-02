using Plots

function boundary_reflection(t_inc, t_bnd)
	MyAngle(tng) = atan(tng[2]/tng[1])
	θ_inc = MyAngle(t_inc)
	θ_bnd = MyAngle(t_bnd)

	c = cos(θ_inc)/t_inc[1]

	θ_inc_flat = θ_inc - θ_bnd
	θ_rfl_flat = -θ_inc_flat
	θ_rfl = θ_rfl_flat + θ_bnd

	return [cos(θ_rfl), sin(θ_rfl)]/c
end

θ_inc = deg2rad(60)
θ_bnd = deg2rad(-45)
t_inc = [cos(θ_inc), sin(θ_inc)]
t_bnd = [cos(θ_bnd), sin(θ_bnd)]
t_rfl = boundary_reflection(t_inc, t_bnd)

x_inc = [-t_inc[1], 0]
y_inc = [-t_inc[2], 0]
x_bnd = [-t_bnd[1], t_bnd[1]]
y_bnd = [-t_bnd[2], t_bnd[2]]
x_rfl = [0, t_rfl[1]]
y_rfl = [0, t_rfl[2]]

pt = plot(aspect_ratio = 1,
	xlim = (-1, 1),
	ylim = (-1, 1))
plot!(x_inc, y_inc)
plot!(x_bnd, y_bnd)
plot!(x_rfl, y_rfl)
display(pt)

θ_rfl = atan(t_rfl[2]/t_rfl[1])

@show rad2deg(θ_rfl + θ_inc)
@show rad2deg(2θ_bnd)