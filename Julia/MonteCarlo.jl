
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
    sample=[]
    while length(sample)<=Nₛ
        a = rand(Uniform(0,ymax))
        b = rand(Uniform(xmin,xmax))
        if a<=pdf(b) push!(sample,b) end
    end
    return sample
end

end

