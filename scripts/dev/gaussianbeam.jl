using OceanAcoustics
scn = ExampleScenarios.n2linear(1)
trc = Trace(scn)
fld = Field(trc)
grid = Grid(fld)
plot_oac(grid)
