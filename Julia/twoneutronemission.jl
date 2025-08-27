# include("/Users/mmadurga/.julia/dev/Utilities/Julia/twoneutronemission.jl

module twoneutronemission

include("/Users/mmadurga/.julia/dev/BetaDecayUtils/src/BetaDecayUtils.jl")

using QuadGK

export gamma3

function gamma3(q,eᵣ₁,eᵣ₂,l1,l2,A)
    lorentz(ε,q,eᵣ,l,A)=BetaDecayUtils.nPenetrability(ε*q,[A-1,1],l)/((ε*q-eᵣ)^2+BetaDecayUtils.nPenetrability(ε*q,[A-1,1],l)^2/4)
    integral, error = quadgk(ε->lorentz(ε,q,eᵣ₁,l1,A)*lorentz(1-ε,q,eᵣ₂,l2,A),0,1)
    return  q*(q-eᵣ₁-eᵣ₂)^2/2π*integral
end

end