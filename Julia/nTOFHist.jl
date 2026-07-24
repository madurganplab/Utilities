module nTOFHist

export εᵥ,
        ToF,
        plotSim,
        MCHist,
        plotMCHist


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

# Rejection sampling is conditional on [t₁,t₂], so its event count must
# be the integral of the plotted rate over that same interval.  A fine
# trapezoidal grid is sufficient for the smooth detector response.
function windowintegral(f,xmin,xmax)
    grid = range(xmin,xmax,length=max(10_001,ceil(Int,(xmax-xmin)*100)+1))
    values = f.(grid)
    step(grid)*(sum(values)-(first(values)+last(values))/2)
end

peakevents = round(Int,ionsample*windowintegral(response,t₁,t₂))
promptevents = round(Int,ionsample*windowintegral(promptresponse,t₁,t₂))
backgroundevents = round(Int,backgroundlevel*(t₂-t₁))

peaksample=MonteCarlo.mcreject(cutoff,t₁,t₂,peakevents,response)

promptsample=MonteCarlo.mcreject(1,t₁,t₂,promptevents,promptresponse)

backgroundsample=MonteCarlo.mcreject(1,t₁,t₂,backgroundevents,t->1.0)

return [ peaksample ; promptsample ; backgroundsample ]

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

end
