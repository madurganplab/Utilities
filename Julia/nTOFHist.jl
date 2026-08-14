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
include("/Users/mmadurga/.julia/dev/Utilities/Julia/FDSiEfficiencies.jl")
include("/Users/mmadurga/.julia/dev/BetaDecayUtils/src/BetaDecayUtils.jl")

using Plots,Measures,Distributions

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

function sampleexcitedbranches(branches,path,intensities,ionsample,t₁,t₂)
    t₁ < t₂ || throw(ArgumentError("t₁ must be smaller than t₂"))
    ionsample >= 0 || throw(ArgumentError("ionsample must be nonnegative"))
    grid = collect(range(t₁,t₂,length=max(10_001,ceil(Int,(t₂-t₁)*100)+1)))
    Δt = step(range(t₁,t₂,length=length(grid)))
    transitionToF = neutronToF.(path,getproperty.(branches,:energy))
    profiles = Vector{Vector{Float64}}(undef,length(branches))
    cdfs = Vector{Vector{Float64}}(undef,length(branches))
    weights = zeros(Float64,length(branches))

    for k in eachindex(branches)
        unitprofile = Float64.(VandleResponses.isolderesponse.(
            grid,1.0,transitionToF[k],path))
        profile = intensities[k].*unitprofile
        profiles[k] = profile
        cumulative = zeros(Float64,length(grid))
        for j in 2:length(grid)
            cumulative[j] = cumulative[j-1] +
                Δt*(unitprofile[j-1]+unitprofile[j])/2
        end
        area = last(cumulative)
        weights[k] = intensities[k]*area
        cdfs[k] = area > 0 ? cumulative./area : cumulative
    end

    eventcount = round(Int,ionsample*sum(weights))
    times = Vector{Float64}(undef,eventcount)
    branchindices = Vector{Int}(undef,eventcount)
    finalindices = Vector{Int}(undef,eventcount)
    if eventcount > 0
        totalweight = sum(weights)
        totalweight > 0 || throw(ArgumentError("total neutron response in [t₁,t₂] is zero"))
        branchcdf = cumsum(weights)./totalweight
        for n in eachindex(times)
            k = searchsortedfirst(branchcdf,rand())
            cdf = cdfs[k]
            u = rand()
            j = searchsortedfirst(cdf,u)
            if j <= 1
                times[n] = first(grid)
            else
                denominator = cdf[j]-cdf[j-1]
                fraction = denominator > 0 ? (u-cdf[j-1])/denominator : 1.0
                times[n] = grid[j-1] + fraction*(grid[j]-grid[j-1])
            end
            branchindices[n] = k
            finalindices[n] = branches[k].final
        end
    end

    return (times=times,branchindices=branchindices,finalindices=finalindices,
            grid=grid,profiles=profiles,transitionToF=transitionToF)
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
            energy > 100eps(Float64) * max(abs(Eₓ[i]), abs(Sₙ), abs(Eᶠ[f]), 1.0) || continue
            maxL = max(8,ceil(Int,Jᵢ[i]+Jᶠ[f]+0.5)+1)
            L = lowestneutronL(Jᵢ[i],πᵢ[i],Jᶠ[f],πᶠ[f]; maxL=maxL)
            width = L > 7 ? 0.0 : BetaDecayUtils.Gamma(energy,A-1,1,L)
            isfinite(width) && width >= 0 ||
                throw(ArgumentError("invalid neutron partial width for transition $i → $f"))
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

function MCHist(cutoff,path,Z,Qᵦ,Sₙ,Eₓ,BGT,backgroundlevel,ionsample,t₁,t₂,eff=nothing;
                returnmetadata=false)

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

sample = [peaksample; promptsample; backgroundsample]
return returnmetadata ?
    (sample=sample,neutrons=peaksample,prompt=promptsample,
     background=backgroundsample) : sample

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
function MCHist(cutoff,path,Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT,
                backgroundlevel,ionsample,t₁,t₂,eff=nothing;
                returnmetadata=false)
    _ = cutoff # Retained for API compatibility; inverse-CDF sampling needs no envelope.
    backgroundlevel >= 0 || throw(ArgumentError("backgroundlevel must be nonnegative"))
    branches = neutronbranches(Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT)
    neutronenergies = getproperty.(branches,:energy)
    intensities = getproperty.(branches,:intensity)
    if isnothing(eff)
        intensities = intensities .* εᵥ.(neutronenergies,0.8)
    end

    neutrondata = sampleexcitedbranches(branches,path,intensities,ionsample,t₁,t₂)
    promptdistribution = Normal(path/c,1)
    promptarea = sqrt(2π)*(cdf(promptdistribution,t₂)-cdf(promptdistribution,t₁))
    promptevents = round(Int,ionsample*promptarea)
    backgroundevents = round(Int,backgroundlevel*(t₂-t₁))
    promptsample = promptevents == 0 ? Float64[] :
        rand(truncated(promptdistribution,t₁,t₂),promptevents)
    backgroundsample = backgroundevents == 0 ? Float64[] :
        rand(Uniform(t₁,t₂),backgroundevents)
    sample = [neutrondata.times; promptsample; backgroundsample]

    return returnmetadata ?
        (sample=sample,neutrons=neutrondata,prompt=promptsample,
         background=backgroundsample,branches=branches,intensities=intensities) : sample
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

simulation = MCHist(cutoff,path,Z,Qᵦ,Sₙ,Eₓ,BGT,backgroundlevel,ionsample,t₁,t₂,eff;
                    returnmetadata=true)
sample = simulation.sample
neutroncount = length(simulation.neutrons)
nbins = max(1,round(Int,t₂-t₁))
binedges = range(t₁,t₂,length=nbins+1)
binwidth = step(binedges)
p=plot(sample,
        seriestype=:stephist,lc=:navy,lw=1.5,bins=binedges,
        weights=fill(inv(binwidth),length(sample)),
        ylims=(0,ymax),xlims=(t₁,t₂),label="",
        xlabel="ToF (ns)",ylabel="counts/ns",
        title="Nₙ = $neutroncount counts",
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

# plot neutron energy scale on top of figure frame
enticks = [0.05,0.1,0.2,0.3,0.5,1.0,2.0,5.0]
tofticks = ToF(path,Qᵦ,Sₙ,enticks.+Sₙ)
for i in findall(tofticks.<(t₂-20))
    plot!([tofticks[i],tofticks[i]],[ymax-ymax/20,ymax],lc=:black,lw=1,label="")
    annotate!(tofticks[i],ymax+ymax/20,text("$(enticks[i])",12))
end
annotate!(t₂,ymax+ymax/20,text("Eₙ (MeV)",12,:right))
return p
end

"""
    plotMCHist(cutoff, ymax, path, Z, A, Qᵦ, Sₙ,
               Eₓ, Jᵢ, πᵢ, Eᶠ, Jᶠ, πᶠ, BGT,
               backgroundlevel, ionsample, t₁, t₂, eff=nothing)

Plot a time-of-flight sample from the excited-state [`MCHist`](@ref) method, its summed expected
response, and every energetically open neutron branch. Parities must be given
as `+1` or `-1`. Returns `(p, stateplots, stateplotsgamma)`, where `p` is the
total spectrum, `stateplots[f]` contains every transition feeding `Eᶠ[f]`
without gamma-efficiency correction, and `stateplotsgamma[f]` contains the
same spectrum corrected by the CLARION gamma efficiency.
"""
function plotMCHist(cutoff,ymax,path,Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT,
                      backgroundlevel,ionsample,t₁,t₂,eff=nothing)
    simulation = MCHist(cutoff,path,Z,A,Qᵦ,Sₙ,Eₓ,Jᵢ,πᵢ,Eᶠ,Jᶠ,πᶠ,BGT,
                        backgroundlevel,ionsample,t₁,t₂,eff;
                        returnmetadata=true)
    sample = simulation.sample
    neutroncount = length(simulation.neutrons.times)
    branches = simulation.branches
    neutrondata = simulation.neutrons
    responsegrid = neutrondata.grid
    profiles = neutrondata.profiles
    totalprofile = reduce(+,profiles; init=zeros(Float64,length(responsegrid)))
    nbins = max(1,round(Int,t₂-t₁))
    binedges = range(t₁,t₂,length=nbins+1)
    binwidth = step(binedges)
    p = plot(sample,
             seriestype=:stephist,lc=:navy,lw=1.5,bins=binedges,
             weights=fill(inv(binwidth),length(sample)),
             ylims=(0,ymax),xlims=(t₁,t₂),label="",
             xlabel="ToF (ns)",ylabel="counts/ns",
             title="Nₙ = $neutroncount counts",
             top_margin=7mm)

    p = plot!(responsegrid,backgroundlevel .+ ionsample.*totalprofile,
              lw=1.5,lc=:red,label="sum")
    p = hline!([backgroundlevel],lw=2,lc=:black,ls=:dash,
               label="background + γ flash")
    p = plot!(t->ionsample*exp(-0.5(t-path/c)^2),0,10,
              lc=:black,lw=2,ls=:dash,label="")

    branchplotindices = collect(1:20:length(responsegrid))
    last(branchplotindices) == length(responsegrid) || push!(branchplotindices,length(responsegrid))
    function plotbranches!(branchplot,branchindices,scale,label)
        isempty(branchindices) && return branchplot
        branchx = Float64[]
        branchy = Float64[]
        sizehint!(branchx,length(branchindices)*(length(branchplotindices)+1))
        sizehint!(branchy,length(branchindices)*(length(branchplotindices)+1))
        for i in branchindices
            append!(branchx,@view responsegrid[branchplotindices])
            append!(branchy,scale.*(@view profiles[i][branchplotindices]))
            push!(branchx,NaN)
            push!(branchy,NaN)
        end
        return plot!(branchplot,branchx,branchy,
                     linecolor=:black,label=label)
    end
    p = plotbranches!(p,eachindex(branches),ionsample,"individual transitions")

    function finalstatespectrum(f,statesample,scale,titlesuffix)
        branchindices = findall(branch->branch.final == f,branches)
        stateprofile = reduce(+,profiles[branchindices];
                              init=zeros(Float64,length(responsegrid)))
        stateplot = plot(statesample,
                         seriestype=:stephist,lc=:navy,lw=1.5,bins=binedges,
                         weights=fill(inv(binwidth),length(statesample)),
                         xlims=(t₁,t₂),label="",
                         xlabel="ToF (ns)",ylabel="counts/ns",
                         title="Eᶠ = $(Eᶠ[f]) MeV$titlesuffix, N = $(length(statesample)) counts",
                         top_margin=7mm)
        stateplot = plot!(stateplot,responsegrid,ionsample*scale.*stateprofile,
                          lw=1.5,lc=:red,label="sum")
        stateplot = plotbranches!(stateplot,branchindices,ionsample*scale,"")
        return stateplot
    end

    stateplots = Plots.Plot[]
    stateplotsgamma = Plots.Plot[]
    for f in eachindex(Eᶠ)
        gammaefficiency = iszero(Eᶠ[f]) ? 1.0 : 2*FDSiEfficiencies.clarion(Eᶠ[f])
        0 <= gammaefficiency <= 1 ||
            throw(ArgumentError("gamma efficiency for Eᶠ=$(Eᶠ[f]) MeV must be between 0 and 1"))
        statesample = neutrondata.times[neutrondata.finalindices .== f]
        gammasample = gammaefficiency == 1 ? copy(statesample) :
            statesample[rand(length(statesample)) .< gammaefficiency]
        push!(stateplots,finalstatespectrum(f,statesample,1.0,""))
        push!(stateplotsgamma,finalstatespectrum(
            f,gammasample,gammaefficiency," (γ efficiency)"))
    end

    # plot neutron energy scale on top of figure frame
    enticks = [0.05,0.1,0.2,0.3,0.5,1.0,2.0,5.0]
    tofticks = neutronToF.(path,enticks)
    function plotneutronenergyscale!(target)
        plotymax = last(ylims(target))
        for i in findall(tofticks.<(t₂-20))
            plot!(target,[tofticks[i],tofticks[i]],
                  [plotymax-plotymax/20,plotymax],
                  lc=:black,lw=1,label="")
            annotate!(target,tofticks[i],plotymax+plotymax/20,
                      text("$(enticks[i])",12))
        end
        annotate!(target,t₂,plotymax+plotymax/20,text("Eₙ (MeV)",12,:right))
        return target
    end
    for target in Iterators.flatten(((p,),stateplots,stateplotsgamma))
        plotneutronenergyscale!(target)
    end

    return p,stateplots,stateplotsgamma
end

end
