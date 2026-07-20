module ShellModelMasses
 
# Equations 8, 9, and 10 from section B of S. Yoshida et al., Phys. Rev. C 97, 054321. 
# DOI: https://doi.org/10.1103/PhysRevC.97.054321

export calculateSMbinding,calculateSMmassExcess,calculateSMmass

"""
calculateSMbinding(Egs,Ecore,Z,N)

Calculate the total shell model binding energy in MeV

Egs:   ground state energy from SM code (MeV, negative value)

Ecore: SM's core binding energy (MeV/A*A, AME2020, positive value)

Z: Atomic number

N: Neutron number
"""
function calculateSMbinding(Egs,Ecore,Z,N)
    A=Z+N
    T=abs(Z-N)/2
    Rc=exp(1.5/A)*A^(1/3)*(0.946-0.573*(2*T/A)^2)
    Ec=0.7*((Z*(Z-1)-0.76*(Z*(Z-1))^(2/3))/Rc)          #Coulomb correction to SM energy
    return -1*(Egs+Ec)+Ecore
end

"""
calculateSMmassExcess(Egs,Ecore,Z,N)

Calculate the shell model mass excess in MeV

Egs:   ground state energy from SM code (MeV, negative value)

Ecore: SM's core binding energy (MeV/A*A, AME2020, positive value)

Z: Atomic number

N: Neutron number
"""
function calculateSMmassExcess(Egs,Ecore,Z,N)
    Δmp,Δmn=7.2890,8.0713                               #proton and neutron mass excess
    bindingenergy=calculateSMbinding(Egs,Ecore,Z,N)     #valence space energy plus core binding energy
    return Z*Δmp+N*Δmn-bindingenergy
end

"""
calculateSMbinding(Egs,Ecore,Z,N)

Calculate the shell model mass in MeV and atomic mass units (u)

Egs:   ground state energy from SM code (MeV, negative value)

Ecore: SM's core binding energy (MeV/A*A, AME2020, positive value)

Z: Atomic number

N: Neutron number
"""
function calculateSMmass(Egs,Ecore,Z,N)
    u,A,Δm=931.49410372,Z+N,calculateSMmassExcess(Egs,Ecore,Z,N)
    return (A*u+Δm),A+Δm/u                              #returns atomic mass in MeV, u
end

end