function linear_interp_fcn(x, y)
	y_interp = linear_interpolation(x, y)
	y_fcn(x) = y_interp(x)
end