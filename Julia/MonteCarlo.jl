
# include("/Users/.julia/dev/Utilities/MonteCarlo.jl")

module MonteCarlo

using Distributions

export mcreject

"""

mcreject(ymax,xmin,xmax,Nₛ,pdf::function)

mcreject implements the Monte Carlo rejection method
(https://web.tecnico.ulisboa.pt/~mcasquilho/acad/theo/simul/Vujic.pdf) 
to create a sample of Nₛ elements from an arbitrary univariate PDF (`pdf`) 
in given `xmin->xmax+xmin` and `0->ymax` ranges.
Requires package "Distributions"

"""
function mcreject(ymax,xmin,xmax,Nₛ,pdf::Function)
    xmin < xmax || throw(ArgumentError("xmin must be smaller than xmax"))
    ymax > 0 || throw(ArgumentError("ymax must be positive"))
    nsamples = round(Int,Nₛ)
    nsamples >= 0 || throw(ArgumentError("Nₛ must be nonnegative"))

    # Check the rejection envelope on a dense grid before sampling.  Without
    # this guard, values above ymax are silently clipped and bias the sample.
    grid = range(xmin,xmax,length=10_001)
    pdfminimum,pdfmaximum = extrema(pdf,grid)
    pdfminimum >= 0 || throw(ArgumentError("pdf must be nonnegative"))
    pdfmaximum <= ymax || throw(ArgumentError(
        "ymax=$ymax is below the sampled PDF maximum $pdfmaximum; " *
        "increase the rejection cutoff"
    ))

    sample=Float64[]
    sizehint!(sample,nsamples)
    while length(sample)<nsamples
        a = rand(Uniform(0,ymax))
        b = rand(Uniform(xmin,xmax))
        pdfvalue = pdf(b)
        pdfvalue <= ymax || throw(ArgumentError(
            "ymax=$ymax is below PDF value $pdfvalue at x=$b; " *
            "increase the rejection cutoff"
        ))
        pdfvalue >= 0 || throw(ArgumentError("pdf must be nonnegative"))
        if a<=pdfvalue push!(sample,b) end
    end
    return sample
end


end
