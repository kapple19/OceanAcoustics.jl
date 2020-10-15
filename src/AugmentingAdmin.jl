function print_properties_types(io::IO, ac::OceanAcoustic)
	println(string(typeof(ac)), "(")
	for p ∈ propertynames(ac)
		println(" ", p, "::", typeof(getproperty(ac, p)))
	end
	print(")")
end