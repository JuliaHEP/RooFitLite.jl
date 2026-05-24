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

function minuitkwargs end
function fitTo end
function visualize end

end
