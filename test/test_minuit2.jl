import Minuit2

@testset "Minuit2 extension" begin
    x = RealVar{Float64}(:x, limits=(0, 10), nbins=50)
    μ = RealVar(:μ, 3.0, limits=(0.0, 5.0))
    σ = RealVar(:σ, 0.8, limits=(0.5, 3.0))
    model = Gaussian(:gauss, x, μ, σ)
    data = generate(model, 2000, nbins=50)

    kwargs = minuitkwargs(model)
    @test haskey(kwargs, :μ)
    @test haskey(kwargs, :limit_μ)

    result = fitTo(model, data)
    @test result isa FitResult
    @test result.engine isa Minuit2.Minuit
    @test result.engine.is_valid
    @test visualize(result.engine, model) isa Plots.Plot
end
