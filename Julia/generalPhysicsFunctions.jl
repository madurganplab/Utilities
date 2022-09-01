# General physics functions and utilities.
# usage --> include("~/generalPhysicsFunctions.jl")

## Spectrometry utilities

#Calculates the gaussian broadening of a 1D spectrum contained in array vectors x coordinate and w Weight)

function gaussianSmoothing(x::Vector,w::Vector,σ)
    smooth=Array{Float64,1}(undef,length(x))
    for i in 1:length(x)
         smooth = smooth .+ w[i]*pdf.(Normal.(x[i],σ),x)
    end
    return smooth
end


## β decay utilities


function daughterActivity(x,A,λ::Float64)
    return A*exp(-λ*x)  #A is initial activity, λ is the decay probability (ln2/T12)
end

function grandDaughterActivity(x,A,λ::Float64,μ::Float64)
    return (λ*μ)/(μ-λ)*A*(exp(-λ*x)-exp(μ*x)) #A:initial activity, λ:daughter, μ:granddaughter
end

function calculateT12(z,Qβ,Ex::Vector,BGT::Vector)
    #Qβ: β decay Q value in MeV
    #Ex: daughter states relative to the ground state energy in MeV
    #z: mother Z 

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
    βEp = (Qβ .-  (Ex)) .* 1000 #convert to keV
    lf = evalCoeff[1] .+ evalCoeff[2].*log.(βEp[findall(βEp.>0)]) .+ evalCoeff[3].*log.(βEp[findall(βEp.>0)]).^2. .+ evalCoeff[4].*log.(βEp[findall(x->x>0,βEp)]).^3.

    D=6144/(-1.2701)^2
    λ=log(2) .* 10 .^lf .* 0.6^2 .* BGT[findall(βEp.>0)] ./ D

    return log(2)./sum(λ)
    
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
    
    sigma = k*R

    s0 = 0
    p0 = sigma
    
    s1 = -1/(1+sigma^2)
    p1 = sigma^3/(1+sigma^2)
    
    s2 = (sigma^2*(2-s1) / ((2-s1)^2+p1^2) ) - 2 
    p2 = p1*(sigma^2/((2-s1)^2+p1^2))

    s3 = (sigma^2*(3-s2) / ((3-s2)^2+p2^2) ) - 3 
    p3 = p2*(sigma^2/((3-s2)^2+p2^2))

    s4 = (sigma^2*(4-s3) / ((4-s3)^2+p3^2) ) - 4
    p4 = p3*(sigma^2/((4-s3)^2+p3^2))

    s5 = (sigma^2*(5-s4) / ((5-s4)^2+p4^2) ) - 5 
    p5 = p4*(sigma^2/((5-s4)^2+p4^2))

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

