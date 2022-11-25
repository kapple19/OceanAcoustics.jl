# Types
abstract type OAC end

function Base.show(io::IO, oac::OAC)
	print(typeof(oac))
	oac isa Depth && print(" (callable)")
	print(":")
	for p in propertynames(oac)
		println()
		prop = getproperty(oac, (p))
		print(" ", p, ": ")
		if typeof(prop) in (Number, String)
			print(prop)
		else
			print(prop |> typeof)
		end
	end
end

# Base.iterate(ocn::Ocean) = ocn

# Base.broadcastable(oac::OAC) = Ref(oac)

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