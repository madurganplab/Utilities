module PhotopeakFit

export gaussianfit,lorentzfit

using LsqFit

"""
gaussianfit(data,xlow,xhigh,param::Vector,n=nothing,lowerbounds=nothing,upperbounds=nothing)

Gamma photopeak fit of a histogram using a gaussian distribution with linear background. 
Returns a three element list containing: parameters from fit; standard error of the parameters; fitting function f(x,param)

data:           histogram in the [energy,counts] format
xlow:           low energy cut for the fit
xhigh:          high energy cut for the fit
param:          initial fit parameters vector (must be 3*n+2). Format: [background-constant,background-slope,area1,centroid1,sigma1,area2,centroid2,sigma2,...]
n:              number of peaks to fit (default=1)
lowerbounds:    lower parameter bounds (optional). same format as param
upperbounds:    upper parameter bounds (optional). same format as param

Example: fit a gamma line at 988 keV

xlow,xhigh,param=980,1000,[600,0.05,1000,988,0.8]
p,s,f = gaussianfit(data,xlow,xhigh,param)

to print the result of the fit parameters

for (i,val) in enumerate(p)
    println("Pi = ",val,"(s[i]))
end

to plot the result of the fit use the returned function f with the optimized parameters p

plot(e->f(e,p),xlow,xhigh)

"""
function gaussianfit(data,xlow,xhigh,param::Vector,n=nothing,lowerbounds=nothing,upperbounds=nothing)
  
    if isnothing(n) n=1 end 

    if length(param)!=n*3+2 throw(DomainError(length(param), "must have [3*number of photopeak + 2] parameters")) end
    if isnothing(lowerbounds) if upperbounds!==nothing throw(DomainError(lowerbounds,"both lower and upper bounds must be defined")) end end
    if isnothing(upperbounds) if lowerbounds!==nothing throw(DomainError(lowerbounds,"both lower and upper bounds must be defined")) end end
    if n>2 throw(DomainError(n, "number of photopeaks must be 2 or less")) end

    fitfunctions = [(e,p)->p[1].+ p[2].*e .+ p[3].*(1 ./sqrt.(2π.*p[5].^2)).*exp.(-0.5.*(e.-p[4]).^2 ./p[5].^2) , 
                    (e,p)->p[1].+ p[2].*e .+ p[3].*(1 ./sqrt.(2π.*p[5].^2)).*exp.(-0.5.*(e.-p[4]).^2 ./p[5].^2) .+ p[6].*(1 ./sqrt.(2π.*p[8].^2)).*exp.(-0.5.*(e.-p[7]).^2 ./p[8].^2)]


    fit=nothing
    
    if isnothing(lowerbounds) && isnothing(upperbounds)
        fit=curve_fit(fitfunctions[n],data[data[:,1].<xhigh .&& data[:,1].>xlow,1],data[data[:,1].<xhigh .&& data[:,1].>xlow,2],param)
    else 
        fit=curve_fit(fitfunctions[n],data[data[:,1].<xhigh .&& data[:,1].>xlow,1],data[data[:,1].<xhigh .&& data[:,1].>xlow,2],param,lower=lowerbounds,upper=upperbounds)
    end

    sigma=stderror(fit)

    return fit.param,sigma,fitfunctions[n]

end


"""
lorentzfit(data,xlow,xhigh,param::Vector,n=nothing,lowerbounds=nothing,upperbounds=nothing)

Gamma photopeak fit of a histogram using a lorentz distribution with linear background. 
Returns a three element list containing: parameters from fit; standard error of the parameters; fitting function f(x,param)

data:           histogram in the [energy,counts] format
xlow:           low energy cut for the fit
xhigh:          high energy cut for the fit
param:          initial fit parameters vector (must be 3*n+2). Format: [background-constant,background-slope,area1,centroid1,sigma1,area2,centroid2,sigma2,...]
n:              number of peaks to fit (default=1)
lowerbounds:    lower parameter bounds (optional). same format as param
upperbounds:    upper parameter bounds (optional). same format as param

Example: fit a gamma line at 988 keV

xlow,xhigh,param=980,1000,[600,0.05,1000,988,0.8]
p,s,f = lorentzfit(data,xlow,xhigh,param)

to print the result of the fit parameters

for (i,val) in enumerate(p)
    println("Pi = ",val,"(s[i]))
end

to plot the result of the fit use the returned function f with the optimized parameters p

plot(e->f(e,p),xlow,xhigh)

"""

function lorentzfit(data,xlow,xhigh,param::Vector,n=nothing,lowerbounds=nothing,upperbounds=nothing)
  
    if isnothing(n) n=1 end 

    if length(param)!=n*3+2 throw(DomainError(length(param), "must have [3*number of photopeak + 2] parameters")) end
    if isnothing(lowerbounds) if upperbounds!==nothing throw(DomainError(lowerbounds,"both lower and upper bounds must be defined")) end end
    if isnothing(upperbounds) if lowerbounds!==nothing throw(DomainError(lowerbounds,"both lower and upper bounds must be defined")) end end
    if n>2 throw(DomainError(n, "number of photopeaks must be 2 or less")) end

    fitfunctions = [(e,p)->p[1].+ p[2].*e .+ p[3]./π.*0.5.*p[5]./((e.-p[4]).^2 .+(0.5.*p[5]).^2) , 
                    (e,p)->p[1].+ p[2].*e .+ p[3]./π.*0.5.*p[5]./((e.-p[4]).^2 .+(0.5.*p[5]).^2) .+ p[6]./π.*0.5.*p[8]./((e.-p[7]).^2 .+(0.5.*p[8]).^2)]


    fit=nothing
    
    if isnothing(lowerbounds) && isnothing(upperbounds)
        fit=curve_fit(fitfunctions[n],data[data[:,1].<xhigh .&& data[:,1].>xlow,1],data[data[:,1].<xhigh .&& data[:,1].>xlow,2],param)
    else 
        fit=curve_fit(fitfunctions[n],data[data[:,1].<xhigh .&& data[:,1].>xlow,1],data[data[:,1].<xhigh .&& data[:,1].>xlow,2],param,lower=lowerbounds,upper=upperbounds)
    end

    sigma=stderror(fit)

    return fit.param,sigma,fitfunctions[n]

end


end