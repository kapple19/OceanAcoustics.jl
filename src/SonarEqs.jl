module SonarEqs

function transmission_loss(pressure::Number)
	TL = -20log10(abs(pressure))
end

end