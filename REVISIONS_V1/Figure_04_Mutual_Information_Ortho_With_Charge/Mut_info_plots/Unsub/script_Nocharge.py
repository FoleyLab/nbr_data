#! /usr/bin/python

import numpy as np
import matplotlib.pylab as pl
import matplotlib as mpl
import sys
import math

if len(sys.argv) != 3:
	print( "Two parameters required: geometry mut_info!")
	exit()

fig, ax = pl.subplots()

# geometry plus pz orbital numbers
geom = np.loadtxt(sys.argv[1],usecols=(1,2,4))

# mutual information
mut_info = np.loadtxt(sys.argv[2])
mi_max = 1.38629436112 

# nalpha
#nalpha = np.loadtxt(sys.argv[3])

# colors
cmap = mpl.cm.cool
#norm = mpl.colors.Normalize(vmin=4, vmax=10)

for i in range(len(geom)):
    for j in range(i + 1, len(geom)):
        d = math.sqrt((geom[i][0] - geom[j][0])*(geom[i][0] - geom[j][0]) + (geom[i][1] - geom[j][1])*(geom[i][1] - geom[j][1]))

		# bond
        if d < 1.9:
            if int(geom[i][2]) == 0 or int(geom[j][2])==0:
                x, y = [geom[i][0], geom[j][0]], [geom[i][1], geom[j][1]]
                pl.plot(x, y, color="k", lw=0.6)
            else:
                x, y = [geom[i][0], geom[j][0]], [geom[i][1], geom[j][1]]
                mi = mut_info[int(geom[i][2])-1][int(geom[j][2])-1]/mi_max
                color = cmap(mi)
                pl.plot(x, y, color=color, lw=3)
                pl.text(np.average(x), np.average(y), f"{mi:.2f}",va="center",ha="center", size = "small",bbox=dict(facecolor='white',edgecolor='none'))

pl.axis('equal')
pl.axis('off')

with open(sys.argv[1]) as fin:
    for line in fin:
        if "#" in line:
            continue
        
        geom = line.split()

        if geom[0]!="C":
            pl.text(float(geom[1]),float(geom[2]),geom[0],va="center",ha="center", size = "xx-large",bbox=dict(facecolor='white',edgecolor='none')) 
#        else:
#            pl.text(float(geom[1]),float(geom[2]),f"{2*nalpha[int(geom[4])-1]:.2f}",va="center",ha="center", size = "small",bbox=dict(boxstyle=f"circle,pad={0.3}",facecolor='white',edgecolor='grey')) 

cax = fig.add_axes([0.1, 0.1, 0.8, 0.05])

#cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
cb1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, orientation='horizontal')
cb1.set_label('Mutual Information/ln(4)')

#pl.axis('off')
#pl.show()
pl.savefig("Plot_Unsubst.pdf")

