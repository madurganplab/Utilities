module DeformationCalculator

using WignerSymbols

export beta2odd,beta2even,qo,cg2

function cg2(Ji::Real, Jf::Real, K::Real)
    # Wigner 3j: (Ji  2  Jf;  K  0  -K)
    threej = wigner3j(Ji, 2, Jf, K, 0, -K)
    phase = (-1)^(Ji - 2 + K)   # phase factor
    return (phase * sqrt(2*Jf + 1) * threej)^2
end

function qo(BE2,CG2)
    return sqrt( 16π/5 * BE2/CG2 )
end

"""

beta2odd(A,Z,Ji,Jf,K,BE2)

β₂ deformation calculator using B(E2) between Ji and Jf members of band K.


"""
function beta2odd(A,Z,Ji,Jf,K,BE2)
    R=1.2*A^(1/3)
    CGsqd=cg2(Ji,Jf,K)
    Q=qo(BE2,CGsqd)
    return sqrt(5π)/3 * Q/(Z*R^2)
end


"""

beta2even(A,Z,BE2↑)

β₂ deformation calculator using B(E2 0⁺→2⁺). For B(E2 2⁺→0⁺) multiply BE2↓ by 5.

"""
function beta2even( A::Int,  Z::Int, BE2::Any) 
    R = 1.2 * A^(1/3)
    return (4π) / (3 * Z * R^2) * sqrt(BE2)
end

end