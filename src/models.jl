abstract type AbstractPdf end

function getproperty(d::AbstractPdf, name::Symbol)
    if hasfield(typeof(d), name)
        return getfield(d, name)
    end

    idx = findfirst(p -> p.name == name, d.params)
    if !isnothing(idx)
        return d.params[idx]
    elseif hasfield(typeof(d), :pdfs)
        idx = findfirst(p -> p.name == name, d.pdfs)
        !isnothing(idx) && return d.pdfs[idx]
    end
    return nothing
end

(model::AbstractPdf)(x) = model.pdf(x, (p.value for p in model.params)...)
show(io::IO, d::AbstractPdf) = print(io, "$(nameof(typeof(d))){$(d.name)} PDF with parameters $([p.name for p in d.params])")

StatsBase.mean(d::AbstractPdf) = StatsBase.mean(distribution(d))
StatsBase.std(d::AbstractPdf) = StatsBase.std(distribution(d))

"""
    generate(d::AbstractPdf, n::Integer=1000; nbins=0)

Generate `n` random numbers from a distribution `d`.
"""
function generate(d::AbstractPdf, n::Int64; nbins=0)
    nbins = nbins > 0 ? nbins : d.x.nbins > 0 ? d.x.nbins : 0
    data = if nbins > 0
        Hist1D(rand(default_rng(), d, n), binedges=range(d.x.limits..., nbins + 1))
    else
        rand(default_rng(), d, n)
    end
    DataSet(data, (d.x,))
end

function Base.rand(rng::AbstractRNG, d::AbstractPdf, n::Int64=1)
    rand(rng, distribution(d), n)
end

"""
    RealVar(name, value=0, error=0; limits=(-Inf, Inf), nbins=0)
    RealVar{T}(name, value=0, error=0; limits=(-Inf, Inf), nbins=0)

Store an observable or floating model parameter with a symbolic `name`, current
`value`, uncertainty `error`, allowed `limits`, and optional histogram `nbins`.
Parameterized construction with `RealVar{T}` creates a constant variable.
"""
mutable struct RealVar{T<:Real}
    const name::Symbol
    value::T
    error::T
    const limits::Tuple{T,T}
    nbins::Int
    const isconst::Bool
end

function RealVar(name, value::T=T(0), error::T=T(0); limits=(-Inf, Inf), nbins=0) where {T<:Real}
    RealVar{T}(name, value, error, limits, nbins, false)
end

function RealVar{T}(name, value=0, error=0; limits=(-Inf, Inf), nbins=0) where {T<:Real}
    RealVar{T}(name, value, error, limits, nbins, true)
end

function setproperty!(v::RealVar{T}, name::Symbol, value) where {T<:Real}
    if name == :value
        v.isconst && throw(ArgumentError("Cannot set value of a constant variable '$(v.name)'"))
        a, b = v.limits
        if value < a
            setfield!(v, name, a)
            @warn("Value $(value) is below limits $(v.limits): limited to $(a)")
        elseif value > b
            setfield!(v, name, b)
            @warn("Value $(value) is above limits $(v.limits): limited to $(b)")
        else
            setfield!(v, name, T(value))
        end
    elseif name == :error
        setfield!(v, name, T(value))
    elseif name == :measurement
        setfield!(v, :value, value(value))
        setfield!(v, :error, uncertainty(value))
    else
        throw(ArgumentError("Field $(name) not found in $(nameof(typeof(v)))"))
    end
end

setindex!(v::RealVar, value) = setproperty!(v, :value, value)
isconst(v::RealVar) = v.isconst

"""
    ConstVar(name, value)
    ConstVar(value=0.0)

Create a constant [`RealVar`](@ref) whose value cannot be changed by a fit.
"""
function ConstVar(name, value::T) where {T<:Real}
    RealVar{T}(name, value, 0, (value, value), 0, true)
end

function ConstVar(value::T=0.0) where {T<:Real}
    RealVar{T}(:none, value, 0, (value, value), 0, true)
end

abstract type AbstractData end

"""
    DataSet(data, observables)

Container for generated or user-provided data together with the observables that
describe the data dimensions.
"""
struct DataSet{T<:Real,N} <: AbstractData
    data::Union{Array{T,N},AbstractHistogram}
    observables::NTuple{N,RealVar{T}}
end

"""
    FitResult(data, model, engine)

Result returned by [`fitTo`](@ref), storing the input data, fitted model, and
underlying minimization engine.
"""
struct FitResult
    data::DataSet
    model::AbstractPdf
    engine
end

"""
    Gaussian(name, x, mean, sigma)

Build a normalized Gaussian PDF in observable `x` with mean and width variables.
The PDF is truncated to `x.limits`.
"""
struct Gaussian{T<:Real,PDF<:Function} <: AbstractPdf
    name::Symbol
    x::RealVar{T}
    μ::RealVar{T}
    σ::RealVar{T}
    params::Tuple{Vararg{RealVar{T}}}
    pdf::PDF
end

function Gaussian(name, x, μ, σ)
    a, b = x.limits
    μ_ = μ.isconst ? μ.value : μ.name
    σ_ = σ.isconst ? σ.value : σ.name
    params = Tuple(p for p in (μ, σ) if !isconst(p))
    fname = gensym("$(name)_pdf")
    pdf = eval(quote
        function $(fname)(x)
            $(fname)(x, $((param.value for param in params)...))
        end
        function $(fname)(x, $((param.name for param in params)...))
            d = _Normal($(μ_), $(σ_))
            scale = cdf(d, $(b)) - cdf(d, $(a))
            ifelse.($(a) .<= x .<= $(b), pdf.(d, x) ./ scale, zero(x))
        end
    end)
    Gaussian(name, x, μ, σ, params, pdf)
end

distribution(d::Gaussian) = truncated(_Normal(d.μ.value, d.σ.value), d.x.limits...)

"""
    Exponential(name, x, c)

Build a normalized exponential PDF in observable `x` with slope parameter `c`.
The PDF is truncated to `x.limits`.
"""
struct Exponential{T<:Real,PDF<:Function} <: AbstractPdf
    name::Symbol
    x::RealVar{T}
    c::RealVar{T}
    params::Tuple{Vararg{RealVar{T}}}
    pdf::PDF
end

function Exponential(name, x, c)
    a, b = x.limits
    c_ = c.isconst ? c.value : c.name
    params = Tuple(p for p in (c,) if !isconst(p))
    fname = gensym("$(name)_pdf")
    pdf = eval(quote
        function $(fname)(x)
            $(fname)(x, $((param.value for param in params)...))
        end
        function $(fname)(x, $((param.name for param in params)...))
            d = _Exponential(-1 / $(c_))
            scale = cdf(d, $(b)) - cdf(d, $(a))
            ifelse.($(a) .<= x .<= $(b), pdf.(d, x) ./ scale, zero(x))
        end
    end)
    Exponential(name, x, c, params, pdf)
end

distribution(d::Exponential) = truncated(_Exponential(-1 / d.c.value), d.x.limits...)

"""
    ArgusPdf(name, x, m0, c, p=ConstVar(:p, 0.5))

Build a normalized ARGUS background PDF in observable `x`.
"""
struct ArgusPdf{T<:Real,PDF<:Function} <: AbstractPdf
    name::Symbol
    x::RealVar{T}
    m₀::RealVar{T}
    c::RealVar{T}
    p::RealVar{T}
    params::Tuple{Vararg{RealVar{T}}}
    pdf::PDF
end

function ArgusPdf(name, m, m₀, c, p=ConstVar(:p, 0.5))
    a, b = m.limits
    m₀_ = m₀.isconst ? m₀.value : m₀.name
    c_ = c.isconst ? c.value : c.name
    p_ = p.isconst ? p.value : p.name
    params = Tuple(p for p in (m₀, c, p) if !isconst(p))
    fname = gensym("$(name)_pdf")
    pdf = eval(quote
        function $(fname)(x)
            $(fname)(x, $((param.value for param in params)...))
        end
        function $(fname)(x, $((param.name for param in params)...))
            d = truncated(ArgusBGDist($(c_), $(p_), 0, $(m₀_)), $(a), $(b))
            ifelse.($(a) .<= x .<= $(b), pdf.(d, x), zero(x))
        end
    end)
    ArgusPdf(name, m, m₀, c, p, params, pdf)
end

distribution(d::ArgusPdf) = truncated(ArgusBGDist(d.c.value, d.p.value, 0, d.m₀.value), d.x.limits...)

"""
    Chebyshev(name, x, coeffs)

Build a Chebyshev PDF in observable `x` from a vector of coefficient variables.
"""
struct Chebyshev{T<:Real,PDF<:Function} <: AbstractPdf
    name::Symbol
    x::RealVar{T}
    coeffs::Vector{RealVar{T}}
    params::Tuple{Vararg{RealVar{T}}}
    pdf::PDF
end

function Chebyshev(name, x, coeffs)
    a, b = x.limits
    params = Tuple(p for p in coeffs if !isconst(p))
    body = "d = ChebyshevDist([" * join([isconst(c) ? "$(c.value)" : "$(c.name)" for c in coeffs], ", ") * "], $(a), $(b))" |> Meta.parse
    fname = gensym("$(name)_pdf")
    pdf = eval(quote
        function $(fname)(x)
            $(fname)(x, $((param.value for param in params)...))
        end
        function $(fname)(x, $((param.name for param in params)...))
            $(body)
            pdf.(d, x)
        end
    end)
    Chebyshev(name, x, collect(coeffs), params, pdf)
end

distribution(d::Chebyshev) = ChebyshevDist([c.value for c in d.coeffs], d.x.limits...)

"""
    AddPdf(name, pdfs, coefs)
    AddPdf(name, pdf1, pdf2, fraction)

Combine PDFs that share one observable. If `coefs` has one entry fewer than
`pdfs`, coefficients are interpreted as recursive fractions. If it has the same
length, coefficients are interpreted as extended yields.
"""
struct AddPdf{T<:Real,PDFS<:Tuple,F<:Function} <: AbstractPdf
    name::Symbol
    x::RealVar{T}
    params::Tuple{Vararg{RealVar{T}}}
    pdfs::PDFS
    fractions::Vector{RealVar{T}}
    pdf::F
    extendable::Bool
    recursive::Bool
end

function AddPdf(name, pdfs, coefs)
    extendable, recursive = false, false
    xs = unique([pdf.x for pdf in pdfs])
    length(xs) != 1 && throw(ArgumentError("All pdfs must have the same x variable $(xs)"))
    params = unique(vcat([param for pdf in pdfs for param in pdf.params], coefs))
    fname = gensym("$(name)_pdf")
    if length(pdfs) == length(coefs)
        extendable = true
        body = join(["$(coef.name) * var\"$(pdf.pdf)\"(x, $(join([param.name for param in pdf.params], ", ")))" for (coef, pdf) in zip(coefs, pdfs)], " + ") |> Meta.parse
    elseif length(pdfs) == length(coefs) + 1
        recursive = true
        coeffs = get_coefficients([fraction.name for fraction in coefs])
        body = join([coeff * " * var\"$(pdf.pdf)\"(x, $(join([param.name for param in pdf.params], ", ")))" for (coeff, pdf) in zip(coeffs, pdfs)], " + ") |> Meta.parse
    else
        throw(ArgumentError("Number of pdfs and coefficients are inconsistent"))
    end
    pdf = eval(quote
        function $(fname)(x, $((param.name for param in params)...))
            $(body)
        end
    end)
    AddPdf(name, xs[1], Tuple(params), Tuple(pdfs), coefs, pdf, extendable, recursive)
end

function AddPdf(name, pdf1::AbstractPdf, pdf2::AbstractPdf, fraction::RealVar)
    AddPdf(name, [pdf1, pdf2], [fraction])
end

show(io::IO, d::AddPdf) = print(io, "$(nameof(typeof(d))){$(join(map(x -> "$(x.name)", d.pdfs), ", "))} PDF with parameters $([p.name for p in d.params])")

function Base.rand(rng::AbstractRNG, d::AddPdf, n::Int64=1)
    r = []
    for i in eachindex(d.pdfs)
        pdf, w = d[i]
        m = Int(round(n * w))
        push!(r, rand(rng, pdf, m))
    end
    return vcat(r...)
end

function get_coefficients(c)
    s = []
    for i in 1:(length(c) + 1)
        f = ""
        for j in 1:(i - 1)
            f *= "(1-$(c[j]))"
            j < length(c) && (f *= " * ")
        end
        push!(s, i > length(c) ? f : f * "$(c[i])")
    end
    return s
end

function getindex(d::AddPdf, idx::Integer)
    len = length(d.pdfs)
    if d.recursive
        f = idx == len ? 1.0 : d.fractions[idx].value
        for i in 1:(idx - 1)
            f *= 1.0 - d.fractions[i].value
        end
        return d.pdfs[idx], f
    else
        total = sum(f.value for f in d.fractions)
        return d.pdfs[idx], d.fractions[idx].value / total
    end
end

function getindex(d::AddPdf, c::Symbol)
    idx = findfirst(pdf -> pdf.name == c, d.pdfs)
    isnothing(idx) && throw(ArgumentError("Component $(c) not found in $(d.name)"))
    return d[idx]
end
