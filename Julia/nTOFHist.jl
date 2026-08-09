module nTOFHist

export εᵥ,
        ToF,
        plotSim,
        MCHist,
        MCHistEx,
        plotMCHist,
        plotMChistEx


include("/Users/mmadurga/.julia/dev/Utilities/Julia/MonteCarlo.jl")
include("/Users/mmadurga/.julia/dev/Utilities/Julia/VandleResponses.jl")
include("/Users/mmadurga/.julia/dev/BetaDecayUtils/src/BetaDecayUtils.jl")

using Plots,Measures

default(framestyle=:box,
        grid=false,
        dpi=600)

c = 30 #cm/ns
mₙ = 998 # MeV/c^2

function εᵥ(Eₙ,εᵦ=0.8)
        calibrationpar = [1.3705E+00, -9.0253E-02, -1.6313E-01, 9.5432E-02, -1.2136E-01, -6.6203E-02, 0]
        if Eₙ<0.05 return 0
        else  return εᵦ*10^sum(calibrationpar.*log10.(Eₙ).^(eachindex(calibrationpar).-1))/100
        end
end

function ToF(path,Qᵦ,Sₙ,Eₓ)
    return  1/sqrt(2) .* path .* sqrt.(mₙ .* (Eₓ[Eₓ.>Sₙ.&&Eₓ.<Qᵦ].-Sₙ)) ./ (c .* (Eₓ[Eₓ.>Sₙ.&&Eₓ.<Qᵦ].-Sₙ))
end

function selecttransitions(values,transitionindices,numberofenergies,name)
    if length(values) == numberofenergies
        return values[transitionindices]
    elseif length(values) == length(transitionindices)
        return values
    else
        throw(ArgumentError("$name must contain either one value per Eₓ energy ($numberofenergies values) or one value per transition between Sₙ and Qᵦ ($(length(transitionindices)) values)"))
    end
end

function branchingfrombgt(z,Qᵦ,Eₓ,BGT)
    betaallowedindices = findall(Eₓ.<Qᵦ)

    if length(BGT) == length(Eₓ)
        betaallowedintensities = BetaDecayUtils.calculateIb(z,Qᵦ,Eₓ,BGT)
    elseif length(BGT) == length(betaallowedindices)
        betaallowedintensities = BetaDecayUtils.calculateIb(z,Qᵦ,Eₓ[betaallowedindices],BGT)
    else
        throw(ArgumentError("BGT must contain either one value per Eₓ energy ($(length(Eₓ)) values) or one value per beta-allowed state below Qᵦ ($(length(betaallowedindices)) values)"))
    end

    Iᵦ = zeros(eltype(betaallowedintensities),length(Eₓ))
    Iᵦ[betaallowedindices] = betaallowedintensities
    return Iᵦ
end

function windowintegral(f,xmin,xmax)
    grid = range(xmin,xmax,length=max(10_001,ceil(Int,(xmax-xmin)*100)+1))
    values = f.(grid)
    return step(grid)*(sum(values)-(first(values)+last(values))/2)
end

neutronToF(path,energy) = path/c*sqrt(mₙ/(2*energy))

function doubledspin(J,name)
    isfinite(J) && J >= 0 || throw(ArgumentError("$name spins must be finite and nonnegative"))
    twoJ = round(Int,2J)
    isapprox(2J,twoJ; atol=1e-10,rtol=0) ||
        throw(ArgumentError("$name spins must be integer or half-integer values"))
    return twoJ
end

function stateparity(parity,name)
    parity in (-1,1) || throw(ArgumentError("$name parities must be +1 or -1"))
    return Int(parity)
end

function lowestneutronL(Ji,πi,Jf,πf; maxL=7)
    twoJi = doubledspin(Ji,"Jᵢ")
    twoJf = doubledspin(Jf,"Jᶠ")
    initialparity = stateparity(πi,"πᵢ")
    finalparity = stateparity(πf,"πᶠ")
    for L in 0:maxL
        initialparity == finalparity*(-1)^L || continue
        # The emitted neutron has j = L ± 1/2 (only j=1/2 for L=0).
        for twoj in unique((abs(2L-1),2L+1))
            lower = abs(twoJf-twoj)
            upper = twoJf+twoj
            lower <= twoJi <= upper && iseven(twoJi-lower) && return L
        end
    end
    throw(ArgumentError("no neutron angular momentum L ≤ $maxL couples Jᵢ=$Ji, πᵢ=$πi to Jᶠ=$Jf, πᶠ=$πf"))
end

function neutronbranches(Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT)
    length(Jᵢ) == length(Eₓ) ||
        throw(ArgumentError("Jᵢ must contain one spin per Eₓ energy ($(length(Eₓ)) values)"))
    length(πᵢ) == length(Eₓ) ||
        throw(ArgumentError("πᵢ must contain one parity per Eₓ energy ($(length(Eₓ)) values)"))
    length(Eᶠ) == length(Jᶠ) ||
        throw(ArgumentError("Eᶠ and Jᶠ must have the same length"))
    length(Eᶠ) == length(πᶠ) ||
        throw(ArgumentError("Eᶠ and πᶠ must have the same length"))
    isempty(Eᶠ) && throw(ArgumentError("Eᶠ must contain at least the N-1 daughter ground state"))
    A > 1 || throw(ArgumentError("A must be greater than 1"))
    all(isfinite,Eᶠ) && all(>=(0),Eᶠ) ||
        throw(ArgumentError("Eᶠ energies must be finite and nonnegative"))
    doubledspin.(Jᵢ,Ref("Jᵢ"))
    doubledspin.(Jᶠ,Ref("Jᶠ"))
    stateparity.(πᵢ,Ref("πᵢ"))
    stateparity.(πᶠ,Ref("πᶠ"))

    Iᵦ = branchingfrombgt(Z,Qᵦ,Eₓ,BGT)
    branches = NamedTuple[]
    for i in eachindex(Eₓ)
        Eₓ[i] < Qᵦ || continue
        openbranches = NamedTuple[]
        for f in eachindex(Eᶠ)
            energy = Eₓ[i]-Sₙ-Eᶠ[f]
            energy > 0 || continue
            maxL = max(8,ceil(Int,Jᵢ[i]+Jᶠ[f]+0.5)+1)
            L = lowestneutronL(Jᵢ[i],πᵢ[i],Jᶠ[f],πᶠ[f]; maxL=maxL)
            width = L > 7 ? 0.0 : BetaDecayUtils.nPenetrability(energy,A-1,1,L)
            isfinite(width) && width >= 0 ||
                throw(ArgumentError("invalid penetrability for transition $i → $f"))
            push!(openbranches,(parent=i,final=f,energy=energy,L=L,width=width))
        end
        isempty(openbranches) && continue
        totalwidth = sum(getproperty.(openbranches,:width))
        totalwidth > 0 || throw(ArgumentError("total neutron width is zero for beta-daughter state $i"))
        for branch in openbranches
            fraction = branch.width/totalwidth
            push!(branches,merge(branch,(fraction=fraction,intensity=Iᵦ[i]*fraction)))
        end
    end
    return branches
end

function plotSim(path,Qᵦ,Sₙ,Eₓ,Iᵦ,t₁,t₂,eff=nothing)
Eₙ = Eₓ[Eₓ.>Sₙ.&&Eₓ.<Qᵦ].-Sₙ
if isnothing(eff) Iₙ = Iᵦ[findall(Eₓ.>Sₙ.&&Eₓ.<Qᵦ)].*εᵥ.(Eₙ,0.8)
else  Iₙ = Iᵦ[findall(Eₓ.>Sₙ.&&Eₓ.<Qᵦ)]
end
    p=plot(t->sum(VandleResponses.isolderesponse.(t,Iₙ,ToF(path,Qᵦ,Sₙ,Eₓ),path)),t₁,t₂,
    lw=2,xlabel="ToF (ns)", ylabel="counts",label="D=$path cm",lc=:red)

for (i,area) in enumerate(Iₙ)
   p=plot!(t->VandleResponses.isolderesponse(t,area,ToF(path,Qᵦ,Sₙ,Eₓ)[i],path),linecolor=:black,label="")
end

display(p)

end

function MCHist(cutoff,path,Z,Qᵦ,Sₙ,Eₓ,BGT,backgroundlevel,ionsample,t₁,t₂,eff=nothing)

peaksample = []
promptsample = []
backgroundsample = []

transitionindices = findall(Eₓ.>Sₙ.&&Eₓ.<Qᵦ)
Eₙ = Eₓ[transitionindices].-Sₙ
Iᵦ = branchingfrombgt(Z,Qᵦ,Eₓ,BGT)
selectedintensities = selecttransitions(Iᵦ,transitionindices,length(Eₓ),"Iᵦ")

if isnothing(eff) Iₙ = selectedintensities.*εᵥ.(Eₙ,0.8)
else  Iₙ = selectedintensities
end

transitionToF = ToF(path,Qᵦ,Sₙ,Eₓ)
response(t) = sum(VandleResponses.isolderesponse.(t,Iₙ,transitionToF,path))
promptresponse(t) = exp(-0.5*(t-path/c)^2)

peakevents = round(Int,ionsample*windowintegral(response,t₁,t₂))
promptevents = round(Int,ionsample*windowintegral(promptresponse,t₁,t₂))
backgroundevents = round(Int,backgroundlevel*(t₂-t₁))

peaksample=MonteCarlo.mcreject(cutoff,t₁,t₂,peakevents,response)

promptsample=MonteCarlo.mcreject(1,t₁,t₂,promptevents,promptresponse)

backgroundsample=MonteCarlo.mcreject(1,t₁,t₂,backgroundevents,t->1.0)

return [ peaksample ; promptsample ; backgroundsample ]

end

"""
    MCHistEx(cutoff, path, Z, A, Qᵦ, Sₙ, Eₓ, Jᵢ, πᵢ, Eᶠ, Jᶠ, πᶠ, BGT,
             backgroundlevel, ionsample, t₁, t₂, eff=nothing)

Generate the same time-of-flight Monte Carlo sample as [`MCHist`](@ref), while
allowing neutron emission to the ground and excited states `Eᶠ` of the
mass-`A-1` daughter. `Jᵢ` and `Jᶠ` are the physical spins (integer or
half-integer) of the beta daughter and neutron daughter states, respectively.
The corresponding parity arrays `πᵢ` and `πᶠ` must contain `+1` or `-1`.

For every energetically open `Jᵢ → Jᶠ` branch, the lowest orbital angular
momentum `L` allowed by angular-momentum coupling is used. Branches from a
given initial state are weighted in proportion to their neutron penetrability,
which corresponds to assuming equal reduced widths. The selected `L` also
satisfies the neutron-emission parity relation `πᵢ = πᶠ(-1)^L`.
"""
function MCHistEx(cutoff,path,Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT,
                  backgroundlevel,ionsample,t₁,t₂,eff=nothing)
    branches = neutronbranches(Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT)
    neutronenergies = getproperty.(branches,:energy)
    intensities = getproperty.(branches,:intensity)
    if isnothing(eff)
        intensities = intensities .* εᵥ.(neutronenergies,0.8)
    end

    transitionToF = neutronToF.(path,neutronenergies)
    response(t) = sum(VandleResponses.isolderesponse.(t,intensities,transitionToF,path))
    promptresponse(t) = exp(-0.5*(t-path/c)^2)

    peakevents = round(Int,ionsample*windowintegral(response,t₁,t₂))
    promptevents = round(Int,ionsample*windowintegral(promptresponse,t₁,t₂))
    backgroundevents = round(Int,backgroundlevel*(t₂-t₁))

    peaksample = MonteCarlo.mcreject(cutoff,t₁,t₂,peakevents,response)
    promptsample = MonteCarlo.mcreject(1,t₁,t₂,promptevents,promptresponse)
    backgroundsample = MonteCarlo.mcreject(1,t₁,t₂,backgroundevents,t->1.0)

    return [peaksample; promptsample; backgroundsample]
end

function plotMCHist(cutoff,ymax,path,Z,Qᵦ,Sₙ,Eₓ,BGT,backgroundlevel,ionsample,t₁,t₂,spin=nothing,eff=nothing)

transitionmask = Eₓ.>Sₙ.&&Eₓ.<Qᵦ
transitionindices = findall(transitionmask)
Eₙ = Eₓ[transitionmask].-Sₙ
Iᵦ = branchingfrombgt(Z,Qᵦ,Eₓ,BGT)
transitionToF = ToF(path,Qᵦ,Sₙ,Eₓ)
selectedintensities = selecttransitions(Iᵦ,transitionindices,length(Eₓ),"Iᵦ")
if isnothing(eff) Iₙ = selectedintensities.*εᵥ.(Eₙ)
else  Iₙ = selectedintensities
end

sample = MCHist(cutoff,path,Z,Qᵦ,Sₙ,Eₓ,BGT,backgroundlevel,ionsample,t₁,t₂,eff)
nbins = max(1,round(Int,t₂-t₁))
binedges = range(t₁,t₂,length=nbins+1)
binwidth = step(binedges)
p=plot(sample,
        seriestype=:stephist,lc=:navy,lw=1.5,bins=binedges,
        weights=fill(inv(binwidth),length(sample)),
        ylims=(0,ymax),xlims=(t₁,t₂),label="",
        xlabel="ToF (ns)",ylabel="counts/ns",
        top_margin=7mm)
        
p=plot!(t->backgroundlevel+sum(VandleResponses.isolderesponse.(t,Iₙ,ToF(path,Qᵦ,Sₙ,Eₓ),path))*ionsample,t₁,t₂,
        lw=1.5,lc=:red,label="sum")
p=hline!([backgroundlevel],lw=2,lc=:black,ls=:dash,label="background + γ flash")
p=plot!(t->ionsample*exp(-0.5(t-path/c)^2),0,10,lc=:black,lw=2,ls=:dash,label="")
if isnothing(spin)
    for (i,area) in enumerate(Iₙ)
        label = i == 1 ? "individual transitions" : ""
        p=plot!(t->VandleResponses.isolderesponse(t,area,transitionToF[i],path)*ionsample,
                linecolor=:black,label=label)
    end
else
    neutronspins = selecttransitions(spin,transitionindices,length(Eₓ),"spin")
    uniquespins = sort(unique(neutronspins))
    categoricalpalette = Plots.palette(:tab10)
    length(uniquespins) <= length(categoricalpalette) ||
        throw(ArgumentError("spin coloring supports at most $(length(categoricalpalette)) unique spins"))
    categoricalcolors = categoricalpalette[1:length(uniquespins)]
    spincolors = Dict(spinvalue => color for (spinvalue,color) in zip(uniquespins,categoricalcolors))
    for spinvalue in uniquespins
        for (transitionnumber,i) in enumerate(findall(==(spinvalue),neutronspins))
            area = Iₙ[i]
            label = transitionnumber == 1 ? "J = $(spinvalue/2)" : ""
            p=plot!(t->VandleResponses.isolderesponse(t,area,transitionToF[i],path)*ionsample,
                    linecolor=spincolors[spinvalue],label=label)
        end
    end
end

enticks = [0.05,0.1,0.2,0.3,0.5,1.0,2.0,5.0]
tofticks = ToF(path,Qᵦ,Sₙ,enticks.+Sₙ)

for i in findall(tofticks.<(t₂-20))
    plot!([tofticks[i],tofticks[i]],[ymax-ymax/20,ymax],lc=:black,lw=1,label="")
    annotate!(tofticks[i],ymax+ymax/20,text("$(enticks[i])",12))
end
annotate!(t₂,ymax+ymax/20,text("Eₙ (MeV)",12,:right))
display(p)

end

"""
    plotMChistEx(cutoff, ymax, path, Z, A, Qᵦ, Sₙ,
                 Eₓ, Jᵢ, πᵢ, Eᶠ, Jᶠ, πᶠ, BGT,
                 backgroundlevel, ionsample, t₁, t₂, eff=nothing)

Plot a time-of-flight sample from [`MCHistEx`](@ref), its summed expected
response, and every energetically open neutron branch. Parities must be given
as `+1` or `-1`.
"""
function plotMChistEx(cutoff,ymax,path,Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT,
                      backgroundlevel,ionsample,t₁,t₂,eff=nothing)
    branches = neutronbranches(Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT)
    neutronenergies = getproperty.(branches,:energy)
    intensities = getproperty.(branches,:intensity)
    if isnothing(eff)
        intensities = intensities .* εᵥ.(neutronenergies,0.8)
    end
    transitionToF = neutronToF.(path,neutronenergies)

    sample = MCHistEx(cutoff,path,Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT,
                      backgroundlevel,ionsample,t₁,t₂,eff)
    nbins = max(1,round(Int,t₂-t₁))
    binedges = range(t₁,t₂,length=nbins+1)
    binwidth = step(binedges)
    p = plot(sample,
             seriestype=:stephist,lc=:navy,lw=1.5,bins=binedges,
             weights=fill(inv(binwidth),length(sample)),
             ylims=(0,ymax),xlims=(t₁,t₂),label="",
             xlabel="ToF (ns)",ylabel="counts/ns",
             top_margin=7mm)

    p = plot!(t->backgroundlevel +
                  sum(VandleResponses.isolderesponse.(t,intensities,transitionToF,path))*ionsample,
              t₁,t₂,lw=1.5,lc=:red,label="sum")
    p = hline!([backgroundlevel],lw=2,lc=:black,ls=:dash,
               label="background + γ flash")
    p = plot!(t->ionsample*exp(-0.5(t-path/c)^2),0,10,
              lc=:black,lw=2,ls=:dash,label="")

    for i in eachindex(branches)
        branch = branches[i]
        label = i == 1 ? "individual transitions" : ""
        p = plot!(t->VandleResponses.isolderesponse(t,intensities[i],transitionToF[i],path)*ionsample,
                  linecolor=:black,label=label)
    end

    enticks = [0.05,0.1,0.2,0.3,0.5,1.0,2.0,5.0]
    tofticks = neutronToF.(path,enticks)
    for i in findall(tofticks.<(t₂-20))
        plot!([tofticks[i],tofticks[i]],[ymax-ymax/20,ymax],lc=:black,lw=1,label="")
        annotate!(tofticks[i],ymax+ymax/20,text("$(enticks[i])",12))
    end
    annotate!(t₂,ymax+ymax/20,text("Eₙ (MeV)",12,:right))
    display(p)
end

end
