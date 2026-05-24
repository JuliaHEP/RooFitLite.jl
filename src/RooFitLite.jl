module RooFitLite

using Distributions
using DistributionsHEP
using FHist
using Random: AbstractRNG, default_rng
using RecipesBase
using StatsBase

import Base: getindex, getproperty, isconst, setindex!, setproperty!, show

const _Exponential = Distributions.Exponential
const _Normal = Distributions.Normal
const _Uniform = Distributions.Uniform
const pdf = Distributions.pdf
const cdf = Distributions.cdf
const truncated = Distributions.truncated
const Hist1D = FHist.Hist1D
const AbstractHistogram = FHist.AbstractHistogram
const wsample = StatsBase.wsample
const ArgusBGDist = DistributionsHEP.ArgusBG
const ChebyshevDist = DistributionsHEP.Chebyshev

export AbstractPdf, AbstractData, RealVar, ConstVar, DataSet, AbstractHistogram, FitResult
export Gaussian, Exponential, ArgusPdf, Chebyshev
export AddPdf, generate, distribution
export minuitkwargs, fitTo, visualize

include("models.jl")
include("recipes.jl")

"""
    minuitkwargs(model; randomize=false)

Return keyword arguments for constructing a Minuit2 minimizer from a RooFitLite
model. This function is provided by the Minuit2 extension.
"""
function minuitkwargs end

"""
    fitTo(model, data)

Fit `model` to `data` and return a [`FitResult`](@ref). This function is
provided by the Minuit2 extension.
"""
function fitTo end

"""
    visualize(args...; kwargs...)

Create plots for models, fit results, or Minuit2 engines when the Plots and
Minuit2 extensions are loaded.
"""
function visualize end

end
