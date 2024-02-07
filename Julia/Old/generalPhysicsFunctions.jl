
# General physics functions and utilities.
# usage --> include("/Users/mmadurga/Library/Mobile Documents/com~apple~CloudDocs/Coding Projects/Julia/generalPhysicsFunctions.jl");

module generalPhysicsFunctions

export logft,calculateT12,gaussianSmoothing,convolution

## Spectrometry utilities

"""

Calculates the gaussian broadening (σ width) of a 1D distribution contained in array vectors μ coordinate and w Weight. 
Projects the distribution to a range defined by xrange
returns vector with distribution normalized to sum(w)

in order to add a background function define:
gaussianSmoothing(x,w,σ,xrange) .+ function.(xrange)

"""
function gaussianSmoothing(μ::Vector,w::Vector,σ,xrange)
    smooth = Array{Float64,1}(undef,length(xrange))
    dx=(xrange[length(xrange)]-xrange[1])/length(xrange)
    for i in 1:length(μ)
         smooth .+= w[i]*pdf.(Normal.(μ[i],σ),xrange).*dx
    end
    return smooth
end

"""

Generic convolution function. 
 Requires a function "f" that takes as variables centroids (from μ) 
 and the evaluation range "xrange" being the x-variable

 example: function test(μ,xrange); return 1 ./ (xrange.-μ).^2; end

"""
function convolution(f::Function,μ::Vector,w::Vector,xrange)
    smooth = Array{Float64,1}(undef,length(xrange))
    dx=(xrange[length(xrange)]-xrange[1])/length(xrange)
    for i in 1:length(μ)
         smooth .+= w[i].* f(μ[i],xrange) .* dx
    end
    return smooth
end

# Generic VANDLE response function (single exponential tail). μ is centroid, σ is gaussian width,
# λ is exponential tail
# in order to use it for the generic convolution routine above it needs to be redefine with fixed σ and λ
# currently NOT normalized to 1

function vandleResponse(μ,σ,λ,xrange)
    vandle = []
    for i in xrange
        if i<μ+5; push!(vandle,exp(-0.5*(i-μ)^2 / σ^2)); end
        if i>=μ+5; push!(vandle,exp(-0.5*(5)^2 / σ^2)*exp(-λ*(i-μ-5))); end
    end  
    return vandle
end

# Advanced VANDLE response function based on Xu's (previously MM's) parameterization.
# Fit is only correct for D=100 cm, with other distances approximated using the fractional method
# sigma'T = (D/D') * sigmaT

function xuresponsefit(t,inputarea,Tₒ,distance)

    Area=Array{Float32,1}(undef,5)

# parameters

    σₒ = (243432. *Tₒ^2. - 11635843. *Tₒ + 736628549.)*1.E-9 
    Aₒ = (0.2*(9961. *Tₒ^2. - 2271025. *Tₒ + 169902047.))*1.E-9
    σ1 = -(32635. *Tₒ^2. + 29611564. *Tₒ - 330764765.)*1E-9
    k1 = (1.3*(3183. *Tₒ^2. - 814553. *Tₒ + 156450993.))*1.E-9
    d1 = (0.65*(71025. *Tₒ^2. - 21045104. *Tₒ + 3969567038.))*1.E-9
    k2 = (1.3*(-1410. *Tₒ^2. + 620420. *Tₒ + 10475713.))*1.E-9
    d2 = (0.7*(79669297. *Tₒ + 11267453221.))*1.E-9
    k3 = (-988. * Tₒ^2. + 315032. * Tₒ - 4530112.)*1.E-9
    d3 = (22333042. * Tₒ^2. - 6632270383. * Tₒ + 539757365929.)*1.E-9

    σEx = (283600. * Tₒ^2. - 1663115. * Tₒ + 681358889.)*1.E-9

    σₒₚ  = sqrt(σₒ^2. +σEx^2.) * 100/distance
    σ₁ₚ = sqrt(σ1^2. +σEx^2.) * 100/distance
    
    b1 = exp(k1*(k1+d1))/(1+(k1+d1)^2.)
    b2 = exp(k1*(k1+d1)-(k1-k2)*(k2+d2))/(1+(k1+d1)^2.)
    b3 = exp(k1*(k1+d1)-(k1-k2)*(k2+d2)-(k2-k3)*(k3+d3))/(1+(k1+d1)^2.)

# conditional variable

    xx = (t-Tₒ)/σ₁ₚ
    
# Calculate Normalization Area 

    Area[1] = -σₒₚ *sqrt(Aₒ)+(1+Aₒ)*σₒₚ *atan(1/sqrt(Aₒ))
    Area[2] = σ₁ₚ * atan(k1+d1)
    Area[3] = σ₁ₚ*(1. -exp(k1*(d1-d2+k1-k2))/(k1*(1. +d1^2. +2. *d1*k1+k1^2.)))
    Area[4] = σ₁ₚ * exp((k1*(d1-d2+k1-k2)-d3*k2-k2*k3))/(k2*(1+d1^2. +2. *d1*k1+k1^2.))*(exp(k2*(d3+k3))-exp(k2*(d2+k2)))
    Area[5] = σ₁ₚ * exp(d1*k1-d2*k1+k1^2+d2*k2-d3*k2-k1*k2+k2^2. -k2*k3)/(k3*(1. +d1^2. +2. *d1*k1+k1^2.))
    
    Norm = sum(Area)
    response = 0.

    # Response Function

    if ((-1/sqrt(Aₒ) <= xx)&&(xx<=0)) 

       response = inputarea * 1. /Norm * (σₒₚ ^2/((t-Tₒ)^2. +σₒₚ ^2.) * (Aₒ+1) - Aₒ )
    
    elseif ((0<xx)&&(xx<=(k1+d1))) 
       
       response = inputarea * 1. /Norm *  σ₁ₚ^2/((t-Tₒ)^2. +σ₁ₚ^2.)
       
    elseif (((k1+d1)<xx)&&(xx<=(k2+d2))) 
       
       response = inputarea * 1. /Norm*b1*exp(-k1*xx)
       
    elseif (((k2+d2)<xx)&&(xx<=(k3+d3))) 
       
       response = inputarea * 1. /Norm*b2*exp(-k2*xx)
       
    elseif (xx>(k3+d3)) then

       response = inputarea * 1. /Norm*b3*exp(-k3*xx)

    end
    
    return response

    end

## β decay utilities


function daughterActivity(x,A,λ::Float64)
    return A*exp(-λ*x)  #A is initial activity, λ is the decay probability (ln2/T12)
end

function grandDaughterActivity(x,A,λ::Float64,μ::Float64)
    return (λ*μ)/(μ-λ)*A*(exp(-λ*x)-exp(μ*x)) #A:initial activity, λ:daughter, μ:granddaughter
end

"""
calculate halflife of the beta decay of an isotope given feedings to excited states

calculateT12(z,Qᵦ,Eₓ,BGT)

Qβ: β decay Q value in MeV

Ex: vector of daughter states relative to the ground state energy in MeV

BGT: vector of BGT values

z: parent Z 

"""
function calculateT12(z,Qᵦ,Eₓ::Vector,BGT::Vector)


    coeff = [ -17.2       7.9015    -2.54        0.28482;
           3.31368   -2.06273    0.703822   -0.075039;
          -0.364018   0.387961  -0.142528    0.016;
           0.0278071 -0.026519   0.0098854  -0.00113772
       ]  ;
    zDaughter = z + 1
    evalCoeff = [
    coeff[1,1] + log(zDaughter) * coeff[1,2] + coeff[1,3]*log(zDaughter)^2. + coeff[1,4]*log(zDaughter)^3.,
    coeff[2,1] + log(zDaughter) * coeff[2,2] + coeff[2,3]*log(zDaughter)^2. + coeff[2,4]*log(zDaughter)^3.,
    coeff[3,1] + log(zDaughter) * coeff[3,2] + coeff[3,3]*log(zDaughter)^2. + coeff[3,4]*log(zDaughter)^3.,
    coeff[4,1] + log(zDaughter) * coeff[4,2] + coeff[4,3]*log(zDaughter)^2. + coeff[4,4]*log(zDaughter)^3.  
    ]
    βEp = (Qᵦ .-  (Eₓ)) .* 1000 #convert to keV
    lf = evalCoeff[1] .+ evalCoeff[2].*log.(βEp[findall(βEp.>0)]) .+ evalCoeff[3].*log.(βEp[findall(βEp.>0)]).^2. .+ evalCoeff[4].*log.(βEp[findall(x->x>0,βEp)]).^3.

    D=6144/(-1.2701)^2
    λ=log(2) .* 10 .^lf .* 0.6^2 .* BGT[findall(βEp.>0)] ./ D

    return log(2)./sum(λ)
    
end

""" 

calculate logft of a given transition to an excitated state

logft(z,t₁₂,Qᵦ,Eₓ,Iᵦ)

Z of the parent, Qᵦ and Eₓ in MeV, t₁₂ in seconds, Iᵦ absolute value

"""
function logft(z,t₁₂,Qᵦ,Eₓ,Iᵦ)  
    coeff = [ -17.2       7.9015    -2.54        0.28482;
    3.31368   -2.06273    0.703822   -0.075039;
   -0.364018   0.387961  -0.142528    0.016;
    0.0278071 -0.026519   0.0098854  -0.00113772
]  ;
zDaughter = z + 1
evalCoeff = [
coeff[1,1] + log(zDaughter) * coeff[1,2] + coeff[1,3]*log(zDaughter)^2. + coeff[1,4]*log(zDaughter)^3.,
coeff[2,1] + log(zDaughter) * coeff[2,2] + coeff[2,3]*log(zDaughter)^2. + coeff[2,4]*log(zDaughter)^3.,
coeff[3,1] + log(zDaughter) * coeff[3,2] + coeff[3,3]*log(zDaughter)^2. + coeff[3,4]*log(zDaughter)^3.,
coeff[4,1] + log(zDaughter) * coeff[4,2] + coeff[4,3]*log(zDaughter)^2. + coeff[4,4]*log(zDaughter)^3.  
]
βEp = (Qᵦ -  (Eₓ)) * 1000 #convert to keV
lf = evalCoeff[1] + evalCoeff[2]*log(βEp) + evalCoeff[3]*log(βEp)^2. + evalCoeff[4]*log(βEp).^3.

return log10(10^lf*t₁₂/Iᵦ)

end

#neutron emission utilities

function nPenetration(x,mass::Vector,Lorb)
    
    AM1 = mass[1] #recoil
    AM2 = mass[2] #neutron
    E = x  * 1000.       #Energy in keV
    R = 1.4  * AM1^0.333333 + AM2^0.333333
    RMAS = AM1*AM2/(AM1+AM2)*931502
    eo = (197329^2)/(2*RMAS*R^2)
    k = sqrt(2*RMAS*E)/197329
    
    σ = k*R

    s0 = 0
    p0 = σ
    
    s1 = -1/(1+σ^2)
    p1 = σ^3/(1+σ^2)
    
    s2 = (σ^2*(2-s1) / ((2-s1)^2+p1^2) ) - 2 
    p2 = p1*(σ^2/((2-s1)^2+p1^2))

    s3 = (σ^2*(3-s2) / ((3-s2)^2+p2^2) ) - 3 
    p3 = p2*(σ^2/((3-s2)^2+p2^2))

    s4 = (σ^2*(4-s3) / ((4-s3)^2+p3^2) ) - 4
    p4 = p3*(σ^2/((4-s3)^2+p3^2))

    s5 = (σ^2*(5-s4) / ((5-s4)^2+p4^2) ) - 5 
    p5 = p4*(σ^2/((5-s4)^2+p4^2))

    if (Lorb==0.) 
        return p0
    elseif (Lorb==1.) 
        return p1
    elseif (Lorb==2.) 
        return p2
    elseif (Lorb==3.) 
        return p3
    elseif (Lorb==4.) 
        return p4
    elseif (Lorb==5.)
        return p5
    end

end

# Probability and sampling utilities

## mcreject implements the Monte Carlo rejection method
## (https://web.tecnico.ulisboa.pt/~mcasquilho/acad/theo/simul/Vujic.pdf) 
## to create a sample of Nₛ elements from an arbitrary univariate PDF (`pdf`) 
## in given `xmin->xmax+xmin` and `0->ymax` ranges.
## Requires package "Distributions"

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
