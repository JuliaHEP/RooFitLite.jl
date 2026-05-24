# RooFitLite.jl

![AI contribution](https://img.shields.io/badge/significant_AI_contribution-human_in_charge-orange.svg)
[![Tests](https://github.com/JuliaHEP/RooFitLite.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/JuliaHEP/RooFitLite.jl/actions/workflows/ci.yml)
[![Documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaHEP.github.io/RooFitLite.jl/dev/)

Lightweight RooFit-style model construction for Julia HEP analyses.

```julia
using RooFitLite

x = RealVar(:x, 0.0, limits=(0.0, 10.0), nbins=50)
μ = RealVar(:μ, 3.0, limits=(0.0, 5.0))
σ = RealVar(:σ, 0.8, limits=(0.5, 3.0))
model = Gaussian(:signal, x, μ, σ)
data = generate(model, 1000)
```

Minuit2 integration is provided as a package extension:

```julia
using RooFitLite
using Minuit2

result = fitTo(model, data)
```
