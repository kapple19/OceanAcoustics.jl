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