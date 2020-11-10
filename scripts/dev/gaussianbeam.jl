using OceanAcoustics
scn = ExampleScenarios.n2linear()
trc = Trace(scn)
fld = Field(trc)
grid = Grid(fld)
plot_oac(grid)
