
# include("/Users/mmadurga/Library/Mobile Documents/com~apple~CloudDocs/Coding Projects/Utilities/Julia/SpectUtil.jl");

module SpectUtil

export gaussianSmoothing,convolution,xuresponsefit

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


"""
xuresponsefit(t,inputarea,Tₒ,distance)

Advanced VANDLE response function based on Xu's (previously MM's) parameterization.
Fit is only correct for D=100 cm, with other distances approximated using the fractional method
sigma'T = (D/D') * sigmaT

"""
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


end
