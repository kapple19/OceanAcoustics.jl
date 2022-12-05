function series125(N::Int)
	@assert N ≥ 1 "Input must be a positive integer."
	series = Int[]
	index = 1
	while isempty(series) || series[end] < N
		for n = [1, 2, 5]
			num = n * index
			if num > N
				return series
			else
				push!(series, num)
			end
		end
		index *= 10
	end
	series
end

export series125