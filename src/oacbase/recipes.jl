@userplot PropagationPlot
@recipe function plot(pp::PropagationPlot)
	fld = pp.args[1]

	legend --> :none
	yflip := true

	@series begin
		seriestype := :heatmap
		seriescolor := cgrad(:jet, rev = true)
		fld.r, fld.z, fld.PL'
	end

end