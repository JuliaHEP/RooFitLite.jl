# RooFitLite.jl

RooFitLite.jl is a small RooFit-style modeling layer for Julia HEP analyses. It
represents observables, parameters, PDFs, generated data, fit results, and plots
as separate Julia objects, following the same compositional idea used by ROOT
RooFit while staying close to ordinary Julia package workflows.

The core workflow is:

1. Define one or more [`RealVar`](@ref) observables and parameters.
2. Build a PDF such as [`Gaussian`](@ref), [`Exponential`](@ref), or
   [`ArgusPdf`](@ref).
3. Generate toy data with [`generate`](@ref), or pass your own data.
4. Fit with [`fitTo`](@ref) when the Minuit2 extension is loaded.
5. Inspect or plot the model and fit result.

```@example quickstart
using RooFitLite
using Random: seed!

seed!(12345)

x = RealVar(:x, 0.0, limits=(0.0, 10.0))
mean = RealVar(:mean, 3.0, limits=(0.0, 5.0))
sigma = RealVar(:sigma, 0.8, limits=(0.1, 3.0))

model = Gaussian(:signal, x, mean, sigma)
data = generate(model, 100)

length(data.data)
```

See the tutorials for a fuller walk through of the RooFit-like pieces and a
signal-plus-background model.
