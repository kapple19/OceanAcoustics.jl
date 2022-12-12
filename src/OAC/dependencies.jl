using Interpolations: linear_interpolation, Line
using IntervalArithmetic: AbstractInterval, Interval, (..)
using ForwardDiff: derivative
using OrdinaryDiffEq: ODEProblem, solve, Tsit5, ContinuousCallback, CallbackSet, terminate!
using RecipesBase: RecipesBase, @userplot, @recipe, @series
using PlotUtils: cgrad
using Statistics: mean

export (..)