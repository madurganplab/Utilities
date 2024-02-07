
# include("/Users/mmadurga/Library/Mobile Documents/com~apple~CloudDocs/Coding Projects/Utilities/Julia/betaDecayUtil.jl");

module bDecayUtil

export logftfromib,calculateT12,daughterActivity, grandDaughterActivity,nPenetrability,logftfrombgt,calculateIb,ecoulomb

## β decay utilities

"""

ecoulomb(Z,N)

calculate the coulomb correction for a shell model calculation of isotope (Z,N)

"""
function ecoulomb(Z,N)
    A = Z + N; T=abs(Z-N)/2
    Rc = ℯ^(1.5/A)*A^(1/3)*(0.946-0.573*(2*T/A)^2)
    return 0.7*(Z*(Z-1)-0.76*(Z*(Z-1))^(2/3))/(Rc)
end

"""

daughterActivity(x,A,λ)

A is initial activity, λ is the decay probability (ln2/T12)

"""
function daughterActivity(x,A,λ::Float64)
    return A*exp(-λ*x)  #A is initial activity, λ is the decay probability (ln2/T12)
end


"""

grandDaughterActivity(x,A,λ,μ)

A:initial daughter activity, λ:daughter, μ:granddaughter

"""
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
    λ=log(2) .* 10 .^lf .* BGT[findall(βEp.>0)] ./ D

    return log(2)./sum(λ)
    
end


"""
calculate branching ratios of the beta decay of an isotope given feedings to excited states

calculateIb(z,Qᵦ,Eₓ,BGT)

Qβ: β decay Q value in MeV

Ex: vector of daughter states relative to the ground state energy in MeV

BGT: vector of BGT values

z: parent Z 

"""
function calculateIb(z,Qᵦ,Eₓ::Vector,BGT::Vector)


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
    λ=log(2) .* 10 .^lf .* BGT[findall(βEp.>0)] ./ D

    return λ./sum(λ)
    
end


""" 

calculate logft of a given transition to an excitated state

logftfromib(z,t₁₂,Qᵦ,Eₓ,Iᵦ)

Z of the parent, Qᵦ and Eₓ in MeV, t₁₂ in seconds, Iᵦ absolute value

"""
function logftfromib(z,t₁₂,Qᵦ,Eₓ,Iᵦ)  
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


"""
logftfrombgt(bgt)

calculate the logft for a given BGT (not quenched)


    
"""
function logftfrombgt(bgt)
    
    return log10((6144.2/(1.2761^2*bgt*0.75^2)))

end


"""

calculates the neutron penetrability p(x,Lorb). Given a general reduced width γ one can calculate the partial width as

Γ = p(x,Lorb)

nPenetrability(x,mass::Vector,Lorb)

x is the excitation energy above Sₙ, Lorb is the neutron amngular momentum

mass[1] is the recoil, mass[2] is the neutron mass (1) in amplitude

"""
function nPenetrability(x,mass::Vector,Lorb)
    
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



end