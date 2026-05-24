using Distributions: Exponential as _Exponential, Normal, pdf, truncated
using FHist
using Plots
using QuadGK

@testset "RooFitLite models" verbose = true begin
    @testset "RealVar" begin
        x = RealVar(:x, 0.0, limits=(0.0, 10.0), nbins=20)
        @test x.name == :x
        @test x.value == 0.0
        @test x.limits == (0.0, 10.0)
        @test x.nbins == 20

        x.value = 1.0
        @test x.value == 1.0
        x[] = 2.0
        @test x.value == 2.0
        @test_throws ArgumentError x.limits = (2.0, 10.0)
        @test_warn "below limits" x.value = -1.0
        @test x.value == 0.0
        @test_warn "above limits" x.value = 11.0
        @test x.value == 10.0

        c = ConstVar(1.0)
        @test c.value == 1.0
        @test_throws ArgumentError c.value = 2.0
    end

    @testset "Gaussian" begin
        x = RealVar(:x, 0.0, limits=(0.0, 10.0))
        μ = RealVar(:μ, 3.0, limits=(0.0, 5.0))
        σ = RealVar(:σ, 0.8, limits=(0.5, 3.0))
        g = Gaussian(:gauss, x, μ, σ)
        @test g.pdf(5.0) ≈ pdf(truncated(Normal(3.0, 0.8), 0.0, 10.0), 5.0)
        @test g.pdf(5.0, 0.0, 1.0) ≈ pdf(truncated(Normal(0.0, 1.0), 0.0, 10.0), 5.0)
        @test quadgk(t -> g.pdf(t), x.limits...)[1] ≈ 1.0

        data = generate(g, 1000)
        @test length(data.data) == 1000
        @test minimum(data.data) >= 0.0
        @test maximum(data.data) <= 10.0
        @test show(devnull, g) === nothing
        @test Plots.plot(g) isa Plots.Plot
        @test visualize(g) isa Plots.Plot
    end

    @testset "Exponential and AddPdf" begin
        x = RealVar{Float64}(:x, limits=(0, 10))
        c = RealVar(:c, -0.5, limits=(-0.8, -0.2))
        e = Exponential(:exp, x, c)
        @test e.pdf(5.0) ≈ pdf(truncated(_Exponential(1 / 0.5), 0.0, 10.0), 5.0)

        μ1 = RealVar(:μ1, 3.0, limits=(0.0, 5.0))
        σ1 = RealVar(:σ1, 0.8, limits=(0.5, 3.0))
        sig1 = Gaussian(:sig1, x, μ1, σ1)
        μ2 = RealVar(:μ2, 6.0, limits=(5.0, 10.0))
        σ2 = RealVar(:σ2, 1.0, limits=(0.5, 3.0))
        sig2 = Gaussian(:sig2, x, μ2, σ2)
        f_sig1 = RealVar(:f_sig1, 0.5, limits=(0.0, 1.0))
        sig = AddPdf(:sig, sig1, sig2, f_sig1)
        @test sig.pdf(5.0, μ1.value, σ1.value, μ2.value, σ2.value, 0.5) ≈ 0.5 * sig1.pdf(5.0) + 0.5 * sig2.pdf(5.0)

        f_sig = RealVar(:f_sig, 0.4, limits=(0.0, 1.0))
        model = AddPdf(:model, sig, e, f_sig)
        @test quadgk(t -> model.pdf(t, (p.value for p in model.params)...), x.limits...)[1] ≈ 1.0
        h = generate(model, 1000, nbins=50)
        @test h.data isa FHist.AbstractHistogram
    end
end
