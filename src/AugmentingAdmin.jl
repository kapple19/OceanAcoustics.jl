function print_properties_types(io::IO, oac::OceanAcoustic)
	println(string(typeof(oac)), "(")
	for p ∈ propertynames(oac)
		println(" ", p, "::", typeof(getproperty(oac, p)))
	end
	print(")")
end