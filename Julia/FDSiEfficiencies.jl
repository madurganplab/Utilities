module FDSiEfficiencies

export clarion,vandle,s3hen

function clarion(Eᵧ)

    # par=[-185.002, 85.8909, -13.0504, 1.12005, -0.116108, 0.00142169, 0.00152923, -0.000100725] #2024 
    par = [-381.802, 154.314, -15.7162, 0.307821, -0.134004, 0.0104717, 0.00243488, -0.000207867] #2025 
    
    return max(sum(par.*log(Eᵧ*1000).^(eachindex(par).-1))/100,0)

end

function vandle(Eₙ)

    par = [-0.6280685283745797,
           -0.08836035581875162,
           -0.19657383508033469,
            0.07603675017078329,
           -0.03366046526582096]

    return max(10^sum(par.*log10(Eₙ).^(eachindex(par).-1)),0)
  
end 

function s3hen(Eₙ)

    par = [-0.11819169450252655
            0.009321834140292796
            -0.05763048959782866
            -0.056427795825639844
            -0.015203721072087794]

    return max(10^sum(par.*log10(Eₙ).^(eachindex(par).-1)),0)        

end

function hagrid(Eᵧ)

    par = [-0.3962390537432245
           -0.6601883326364469
            0.06807068383125314
            0.060087348759686306
           -0.4294040973226941
    ]

    Ωhagrid = 20*π*(5.08/2)^2/(4*π*18.8^2); # 2 inch

    return Ωhagrid *  max(10^sum(par.*log10(Eᵧ).^(eachindex(par).-1)),0)

end


end