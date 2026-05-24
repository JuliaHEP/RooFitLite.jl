const plot_attributes = Dict{Symbol,Any}()

info(model) = join(["$(p.name) = $(round(p.value, sigdigits=3)) ± $(round(p.error, sigdigits=3))" for p in model.params], '\n')

@recipe function f(pdf::AbstractPdf; components=())
    seriestype --> :path
    scale = get(plotattributes, :integral, get(plot_attributes, :integral, 1.0))
    @series begin
        pdf_int = pdf isa AddPdf && pdf.extendable ? sum(f.value for f in pdf.fractions) : 1
        linestyle := :solid
        components != () && (label := string(pdf.name))
        (x -> pdf.pdf(x, (p.value for p in pdf.params)...) * scale / pdf_int, pdf.x.limits...)
    end
    for c in components
        comp, weight = pdf[c]
        a, b = comp.x.limits
        func = x -> comp.pdf(x, (p.value for p in comp.params)...) * scale * weight
        @series begin
            label := string(comp.name)
            linestyle := :dash
            (func, a, b)
        end
    end
end

@recipe function f(ds::DataSet)
    seriestype --> :scatter
    x = ds.observables[1]
    label --> string(x.name)
    ylabel --> "events"
    markercolor --> :black
    if ds.data isa AbstractHistogram
        h = ds.data
        bins = h |> nbins
    else
        bins = get(plotattributes, :bins, 0)
        bins = bins > 0 ? bins : x.nbins > 0 ? x.nbins : 100
        h = Hist1D(ds.data, binedges=range(x.limits..., bins + 1))
    end
    x := bincenters(h)
    y := bincounts(h)
    yerr := sqrt.(bincounts(h))
    plot_attributes[:bins] = bins
    plot_attributes[:integral] = integral(h, width=true)
    ()
end

@recipe function f(r::FitResult)
    (; data, model) = r
    legend --> :outerleft
    legendfontsize --> 6
    labelfontsize --> 8
    legend_background_color --> :lightgray
    ylabel --> "events"
    xlabel --> string(r.model.x.name)
    @series begin
        label := "data"
        color := :black
        (data,)
    end
    @series begin
        label --> info(model)
        (model,)
    end
end
