# Functions utilities package for Chi squared minimization of binned distributions
# Functions take a parameter vector of any length. Define the vector of the approprite length
# for Optim to minimize the chi squared function
# All functions require a global `x` vector containing the bin values for the sample being fitted
# usage --> include("PhysicsFitFunctions.jl")

#Classical gaussian (normal distribution). Para[3] corresponds to amplitude, not area.
function gauss(para::Vector) 
    for i in [1:length(x)]; return para[3] .* exp.(-0.5.*(para[1] .- collect(x)[i]).^2 /para[2]^2) end  
 end

 #Exponential decay, requires one decay constant para[1] and one initial sample size para[2]
 function expFit(para::Vector) 
    for i in [1:length(x)]; return para[2] .* exp.(-1 .*para[1].*collect(x)[i]) end  
 end

 #Exponential decay of the daughter, requires two decay parameters, lambda::para[1] and mu::para[2], and one mother initial sample size para[3]
 function daughterExpFit(para::Vector) 
    for i in [1:length(x)]; return (para[1].*para[2])./(para[2].-para[1]).*para[3].*(exp.(-para[1].*collect(x)[i])-exp.(-para[2].*collect(x)[i])) end  
 end

end