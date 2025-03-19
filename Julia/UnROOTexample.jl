using Pkg; Pkg.add("UnROOT")
using UnROOT,Plots

default(lc=:navy,
        lw=1.5,
        xlabel="E(keV)",
        tickfontsize=11,
        framestyle=:box,
        grid=false)

f=ROOTFile("/SCRATCH/DScratch3/is733/2024apr/K47_T1_001_DD.root")

t=LazyTree(f,"PixTree",[r"clover_vec",r"singlebeta"])



#Plots of singles and betagated clover energies


cloverenergy=[]
betagated=[]

for event in t

        betanum=event.PixTreeEvent_singlebeta_vec__singlebeta_vec__detNum
        energy=event.PixTreeEvent_clover_vec__clover_vec__energy

        append!(cloverenergy,energy)

        if length(betanum) > 0
                append!(betagated,energy)
        end

end

histogram(cloverenergy);

savefig("singles.png")

histogram(betagated)

savefig("betagated.png")

histogram(cloverenergy)
histogram!(betagated,lc=:red)

savefig("comparison.png")









