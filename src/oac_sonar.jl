module SonarEqs

function transmission_loss(pressure::Number)
	TL = -20log10(abs(pressure))
end

# Export all
for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
        @eval export $n
    end
end

end