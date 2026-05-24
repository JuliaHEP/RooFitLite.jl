module RooFitLiteMinuit2Ext

using FHist
using Minuit2
using RooFitLite

function RooFitLite.minuitkwargs(d::RooFitLite.AbstractPdf; randomize=false)
    value(p) = if randomize
        mul = rand(RooFitLite._Uniform(-0.1, 0.1))
        p.value + mul * (mul > 0 ? (p.limits[2] - p.value) : (p.value - p.limits[1]))
    else
        p.value
    end
    merge(
        Dict(p.name => value(p) for p in d.params),
        Dict(Symbol(:limit_, p.name) => p.limits for p in d.params),
    )
end

function RooFitLite.fitTo(model::RooFitLite.AbstractPdf, ds)
    data = ds isa RooFitLite.DataSet ? ds.data : ds
    cost = if data isa FHist.AbstractHistogram
        if !isnothing(model.extendable) && model.extendable
            Minuit2.ExtendedBinnedNLL(data, model.pdf, use_pdf=:approximate)
        else
            Minuit2.BinnedNLL(data, model.pdf, use_pdf=:approximate)
        end
    elseif data isa AbstractArray
        if !isnothing(model.extendable) && model.extendable
            fname = gensym("$(model.name)_ext_pdf")
            ext_pdf = eval(quote
                function $(fname)(x, $((param.name for param in model.params)...))
                    +($((f.name for f in model.fractions)...)), $(model.pdf)(x, $((param.name for param in model.params)...))
                end
            end)
            Minuit2.ExtendedUnbinnedNLL(data, ext_pdf)
        else
            Minuit2.UnbinnedNLL(data, model.pdf)
        end
    else
        throw(ArgumentError("Data type not supported"))
    end

    m = Minuit2.Minuit(cost; RooFitLite.minuitkwargs(model)...)
    Minuit2.migrad!(m)
    if m.is_valid
        for (param, value, error) in zip(model.params, m.values, m.errors)
            param.value = value
            param.error = error
        end
        return RooFitLite.FitResult(ds, model, m)
    else
        throw(ArgumentError("Fit failed"))
    end
end

end
