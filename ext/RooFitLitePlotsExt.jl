module RooFitLitePlotsExt

using Plots
using RooFitLite

function RooFitLite.visualize(model::RooFitLite.AbstractPdf; kwargs...)
    isnothing(model.x) && throw(ArgumentError("Model does not have a variable"))
    return Plots.plot(x -> model(x), model.x.limits...; label="$(model.name)", kwargs...)
end

end
