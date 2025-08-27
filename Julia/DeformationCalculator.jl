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

function beta2odd(A,Z,Ji,Jf,K,BE2)
    R=1.2*A^(1/3)
    CG=cg2(Ji,Jf,K)
    Q=qo(BE2,CG)
    return sqrt( 5π)/3 * Q/(Z*R^2 )
end

function beta2even( A::Int,  Z::Int, BE2::Float64)
    R = 1.2 * A^(1/3)
    return (4π) / (3 * Z * R^2) * sqrt(BE2)
end



end