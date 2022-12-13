module OACBase
using PlotUtils: cgrad
using RecipesBase: RecipesBase #= check =#, @userplot, @recipe, @series
using IntervalArithmetic: AbstractInterval, Interval, (..)
using Interpolations: linear_interpolation, Line

include("base/types.jl")
include("base/auxiliaries.jl")
include("base/scenarios.jl")

export (..)
end # module OACBase

using .OACBase

# Re-Exports
export (..)

# Exports
## Exceptions
export NotSorted
export NotAllUnique

## Types
export Surface
export Bottom
export Ocean
export Environment
export Source
export Receiver
export Entities
export Scenario
export Field

## Auxiliary Functions
export calc_ocean_depth_range

## Plot Recipes
export scenarioplot
export scenarioplot!
export propagationplot