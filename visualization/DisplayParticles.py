#!/usr/bin/env python

##
##--- 
##

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, sys, csv

## Parameters visualisation
shouldPlotEnd = False
saveVideo     = True
maxN          = 2000            # max number of particles

##------------------ 0.1) Read parameters ------------------##
##----------------------------------------------------------##
f = open('../bin/PARAMETER_MicroVic_flat.txt');
lines=f.readlines()
f.close()
Lx        = float(lines[12-1])
Ly        = float(lines[13-1])
BC        =   int(lines[16-1])
dt        = float(lines[18-1])
Time      = float(lines[19-1])
nSeed     =   int(lines[27-1])
dx        = float(lines[31-1])
dxy       = float(lines[33-1])
jumpPrint =   int(lines[35-1])

##------------------ 0.2) Initialisation ------------------##
##---------------------------------------------------------##
nTime = int(Time/dt+.5);

def loadBinary(nameFile):
    # Read the binary data in the file 'nameFile' in 'l'.
    # We also need to precise the number of rows (numberRow).
    f = open(nameFile, "rb") 
    f.seek(4, os.SEEK_SET)
    l = np.fromfile(f, dtype=np.float64)
    f.close()
    return l
if (shouldPlotEnd):
    jumpPrint = nTime;

# l = np.fromfile(nameFile, dtype=np.float32, count=-1, sep='')
# l = np.fromfile(nameFile, dtype=float, count=-1, sep='')
    
# init
#-----
if (saveVideo):
    plt.ioff()
else:
    plt.ion()
# plt.show()
plt.clf()


##---------------------------------------------------------------##
##------------------------  la boucle  --------------------------##
##---------------------------------------------------------------##

#- A) first image
# load
#extension = 'seed' + str(nSeed) + '_' + str(0).zfill(9) + '.udat'
extension = str(0).zfill(9) + '.udat'
x     = loadBinary('../data/particleX_'     + extension)[0:maxN]
y     = loadBinary('../data/particleY_'     + extension)[0:maxN]
theta = loadBinary('../data/particleTheta_' + extension)[0:maxN]
# plot
pltParticle  = plt.quiver(x,y,
                          np.cos(theta),np.sin(theta),scale=40,
                          headaxislength=8,color='b',linewidth=.2,pivot='mid') # !!! pivot='mid','tip'
# deco
fig1 = plt.plot()
#plt.axis([0, Lx, 0, Ly])
plt.axis([0-.1, Lx+.1, 0-.1, Ly+.1])
plt.axes().set_aspect('equal', 'box')
plt.xlabel(r'$x$',fontsize=25)
plt.ylabel(r'$y$',fontsize=25)
theTitle = plt.title("Time t=" +'{:04.2f}'.format(0*dt),
                     horizontalalignment='center',fontsize=20)
# extra deco (boundary condition)
if (BC==2):
    # box
    rect = mpatches.Rectangle(
        (0.0, 0.0),   # (x,y)
        Lx,Ly,         # width, height
        alpha=.1,color='r')
    rect2 = mpatches.Rectangle(
        (0.0, 0.0),
        Lx,Ly,
        fill=False,lw=3,color='r')
    plt.gca().add_patch(rect)
    plt.gca().add_patch(rect2)
elif (BC==4):
    # channel
    rect  = mpatches.Rectangle((0.0-0.1,0.0),.1,Ly,color='r')
    rect2 = mpatches.Rectangle((Lx+0.0,0.0),.1,Ly,color='r')
    plt.gca().add_patch(rect)
    plt.gca().add_patch(rect2)
elif (BC==5):
    # circle
    circle  = mpatches.Circle((Lx/2,Ly/2),Lx/2,color='r',fill=None,lw=3)
    circle2  = mpatches.Circle((Lx/2,Ly/2),Lx/2,color='r',alpha=.1,lw=0)
    plt.gca().add_patch(circle)
    plt.gca().add_patch(circle2)
plt.draw()

#- B)  boucle
for iTime in range(0,nTime+1,jumpPrint):
    # load
    #extension = 'seed' + str(nSeed) + '_' + str(iTime).zfill(9) + '.udat'
    extension = str(iTime).zfill(9) + '.udat'
    x     = loadBinary('../data/particleX_'     + extension)[0:maxN]
    y     = loadBinary('../data/particleY_'     + extension)[0:maxN]
    theta = loadBinary('../data/particleTheta_' + extension)[0:maxN]
    # update
    pltParticle.set_offsets(np.c_[x,y])
    pltParticle.set_UVC(np.cos(theta),np.sin(theta))
    # deco
    theTitle.set_text('Time t=' +'{:04.2f}'.format(iTime*dt))
    plt.draw()
    if (saveVideo):
        plt.savefig('images/particles_'+str(iTime).zfill(9) +'.jpg', bbox_inches='tight')

##---------------------------------------------------------------##
##---------------------------------------------------------------##

## mencoder ''mf://*.jpg'' -mf fps=40 -o ../bidon3.avi -ovc lavc -lavcopts vcodec=mpeg4
## mencoder ''mf://*.jpg'' -mf fps=25 -o ../bidon3_bis.avi -ovc x264 -lavcopts vcodec=mpeg4
## mencoder ''mf://*.jpg'' -mf fps=40 -o ../output3.avi -oac copy -ovc copy

## mencoder ''mf://*.jpg'' -mf fps=25 -o ../bidon3_bis.avi -ovc x264 -lavcopts vcodec=mpeg4 -vf scale=640:-2
