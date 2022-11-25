# Types
abstract type OAC end

function Base.show(io::IO, oac::OAC)
	print(typeof(oac))
	oac isa Depth && print(" (callable)")
	println(":")
	for p in propertynames(oac)
		prop = getproperty(oac, (p))
		print(" ", p, ": ")
		if prop isa Number
			println(prop)
		else
			println(prop |> typeof)
		end
	end
end

# Exceptions
struct NotSorted <: Exception
	var
end
Base.showerror(io::IO, e::NotSorted) = print(io, e.var, " not sorted")
export NotSorted

struct NotAllUnique <: Exception
	var
end
Base.showerror(io::IO, e::NotAllUnique) = print(io, e.var, " not all unique")
export NotAllUnique