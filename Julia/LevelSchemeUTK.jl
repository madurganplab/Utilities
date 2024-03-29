# This program (levelSchemeUTK) is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# Support and maintenance NOT GUARANTEED

### Dr. Miguel Madurga, University of Tennessee (2021)

using Pkg
Pkg.add("Plots")


using Plots, LaTeXStrings

#  define annotations for the levels dictionary
global annotationsFonts=[ Plots.font(14,"helvetica",:left,:black),     # old level jpi [1]
                   Plots.font(9,"helvetica",:right,:black),    # old level E [2]
                   Plots.font(14,"helvetica",:left,:red),       # new level jpi [3]
                   Plots.font(9,"helvetica",:right,:red),      # new level Energy [4]
                   Plots.font(9,"helvetica",:right,:black),     # old level lifetime [5]
                   Plots.font(9,"helvetica",:right,:red),        # new level lifetime [6]
                   Plots.font(9,"helvetica oblique",:left,:black), #normal gamma line [7]
                   Plots.font(9,"helvetica oblique",:left,:red),               #new gamma line [8]
                   Plots.font(18,"helvetica",:center,:black)    #Nucleus] [9]
                   ]                 

annotationsFonts[7].rotation=70 ; annotationsFonts[8].rotation=70  #gamma transition label angle

# # Dictionary containing levels and transitions information

# # Energy-1 LineColor-2 LineWidth-3 LabelJpi-4 LabelFont-5 LabelE-6 LabelFont-7 labelLifetime-8 lifetimeFont-9
levels=[ 
0    :black 3 L"0^+"   annotationsFonts[1]   " "        annotationsFonts[2]  " xx ms "           annotationsFonts[5] ; #1
1000  :black 1 L"2^+"   annotationsFonts[1]   "1000 keV"  annotationsFonts[2]  " xx ps "         annotationsFonts[5] ; #2
1500  :red   1 L"4^+" annotationsFonts[3]   "1500 keV"  annotationsFonts[4]  "  "                annotationsFonts[6] ;  #3
2000 :black 1 L"6^+"   annotationsFonts[1]   "2000  keV" annotationsFonts[2]  "  "               annotationsFonts[5]    #4
] 

# # levelHigh-1 levelLow-2 arrowColor-3 arrowWidth-4 Label-5 labelFont-6
transitions=[
levels[2,1] levels[1,1] :black 1 "1000 keV [E2]"  annotationsFonts[7] ;
levels[3,1] levels[2,1] :red   1 "500 keV"       annotationsFonts[8] ;
levels[4,1] levels[3,1] :black 1 "500 keV"      annotationsFonts[7]
]

# Canvas and position settings
canvasWidth = 80;
canvasLow = -170             #in keV 
canvasHigh = 2400            #in keV
goldenRatio = 20/(2400+170); nucLabel = canvasLow + goldenRatio*(canvasHigh-canvasLow)

levelLineStart = canvasWidth*0.25
levelLineStop = canvasWidth*0.75

energyShift=90              #define annotation y position level (keV) 
lineShift=30                #define annotation y position gamma line (keV)
canvasShift=12              #position of initial gamma line with respect of the canvas edge (points)
separation=5               #separation between gamma lines

# # the one time we need to define an empty canvas to define the drawing region...

plot([],ticks=:false,grid=:false,label=:false,axis=:false,
ylims=(canvasLow,canvasHigh),xlims=(00,canvasWidth)
)

# # Nuclei labels
annotate!([(canvasWidth/2,nucLabel,(L"^{A}_{Z}\textrm{Uub}_{N}",annotationsFonts[9]))])

# # Plot levels and transitions from dictionaries

for i  in eachrow(levels)
    plot!([levelLineStart,levelLineStop], [i[1],i[1]],
    lcolor=i[2],label=:false,linewidth=i[3]
    )
    annotate!([(levelLineStart,i[1]+energyShift,(i[4],i[5])),
    (levelLineStop,i[1]+energyShift-30,(i[6],i[7])),
    (levelLineStop+length(i[8]),i[1],(i[8],i[9]))])
end

for j  in eachrow(transitions)

    plot!([levelLineStop-canvasShift,levelLineStop-canvasShift],[j[1],j[2]],arrow=true,label=false,linecolor=j[3],lw=j[4])
    annotate!([ (levelLineStop-canvasShift,j[1]+lineShift,(j[5],j[6])) ])
    global canvasShift += separation
 
end

# last plot needed to force overlay from for-loops (the loop exists with null, that is, no plot!)
plot!()

