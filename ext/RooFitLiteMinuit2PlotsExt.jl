module RooFitLiteMinuit2PlotsExt

using Minuit2
using Plots
using RooFitLite
using FHist

function RooFitLite.visualize(m::Minuit2.Minuit, model::RooFitLite.AbstractPdf; nbins=-1, components=(), kwargs...)
    isnothing(m.cost) && throw(ArgumentError("Minuit object does not have a cost function"))
    m.cost.ndim > 1 && throw(ArgumentError("Cost function dimension > 1 not supported"))
    m.cost isa Minuit2.UnbinnedCostFunction && nbins == -1 && (nbins = 50)
    plt = visualize_cost(m.cost, m.is_valid, collect(m.values), nbins=nbins; kwargs...)
    for c in components
        comp, weight = model[c]
        a, b = comp.x.limits
        scale = if m.cost isa Minuit2.UnbinnedNLL
            nbins == -1 && (nbins = 50)
            weight * prod(Base.size(m.cost.data)) * (b - a) / nbins
        elseif m.cost isa Minuit2.ExtendedUnbinnedNLL
            nbins == -1 && (nbins = 50)
            weight * (b - a) / nbins
        elseif m.cost isa Minuit2.BinnedNLL
            nbins != -1 && nbins != comp.x.nbins && @warn "Forced #bins to $(comp.x.nbins)"
            weight * (b - a) / comp.x.nbins * sum(m.cost.bincounts)
        elseif m.cost isa Minuit2.ExtendedBinnedNLL
            nbins != -1 && nbins != comp.x.nbins && @warn "Forced #bins to $(comp.x.nbins)"
            weight * (b - a) / comp.x.nbins
        end
        func = x -> comp.pdf(x, (p.value for p in comp.params)...) * scale
        Plots.plot!(plt, func; label="$(comp.name)", kwargs...)
    end
    return plt
end

function visualize_cost(cost::Minuit2.LeastSquares, is_valid, pars; nbins=50, kwargs...)
    x = cost.x
    y = cost.y
    yerr = cost.yerror
    plt = Plots.plot(x, y, yerr=yerr, seriestype=:scatter, kwargs...)
    if is_valid
        yt = cost.model.(x, pars...)
        Plots.plot!(plt, x, yt; label="Fit")
    end
    return plt
end

function visualize_cost(cost::Minuit2.UnbinnedCostFunction, is_valid, pars; nbins=50, kwargs...)
    h = FHist.Hist1D(cost.data, nbins=nbins)
    x = FHist.bincenters(h)
    y = FHist.bincounts(h)
    dy = sqrt.(y)
    plt = Plots.plot(x, y; yerr=dy, seriestype=:scatter, label="Data", kwargs...)
    if is_valid
        if cost isa Minuit2.UnbinnedNLL
            scale = prod(Base.size(cost.data)) * (x[2] - x[1])
            Plots.plot!(plt, x -> cost.model(x, pars...) * scale; label="Fit")
        else
            scale = x[2] - x[1]
            Plots.plot!(plt, x -> cost.model(x, pars...)[2] * scale; label="Fit")
        end
    end
    return plt
end

function visualize_cost(cost::Minuit2.BinnedCostFunction, is_valid, pars; nbins=50, kwargs...)
    x = cost.bincenters
    dx = (x[2] - x[1]) / 2
    y = cost.bincounts
    dy = sqrt.(y)
    plt = Plots.plot(x, y; yerr=dy, seriestype=:scatter, label="Data", kwargs...)
    if is_valid
        scale = cost isa Minuit2.BinnedNLL ? sum(cost.bincounts) : 1
        if cost.use_pdf == :approximate
            f = x -> cost.model(x, pars...) * scale * 2dx
        else
            f = x -> (cost.model(x + dx, pars...) - cost.model(x - dx, pars...)) * scale
        end
        Plots.plot!(plt, f; label="Fit")
    end
    return plt
end

function visualize_cost(cost::Minuit2.CostSum, is_valid, pars; nbins=50, kwargs...)
    plots = []
    for (i, c) in enumerate(cost.costs)
        plt = visualize_cost(c, is_valid, pars[cost.argsmapping[i]]; nbins=nbins, kwargs...)
        push!(plots, plt)
    end
    n = plots |> length |> sqrt |> ceil |> Int
    return Plots.plot(plots...; layout=(n, n))
end

end
