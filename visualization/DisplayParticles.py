#!/usr/bin/env python

##
##--- 
##

import time
import numpy as np
import pylab as P
import matplotlib.pyplot as plt


import os, sys, csv

## Parameters visualisation
shouldPlotEnd = False
saveVideo     = False
lengthArrow   = 1;


##------------------ 0.1) Read parameters ------------------##
##----------------------------------------------------------##
f = open('../bin/PARAMETER_MicroVic_flat.txt');
lines=f.readlines()
f.close()
Lx        = float(lines[12-1])
Ly        = float(lines[13-1])
dt        = float(lines[18-1])
Time      = float(lines[19-1])
nSeed     = int(lines[27-1])
dx        = float(lines[31-1])
dxy       = float(lines[33-1])
jumpPrint = int(lines[35-1])

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
x     = loadBinary('../data/particleX_'     + extension)
y     = loadBinary('../data/particleY_'     + extension)
theta = loadBinary('../data/particleTheta_' + extension)
# plot
pltParticle  = plt.quiver(x,y,
                          lengthArrow*np.cos(theta),lengthArrow*np.sin(theta),
                          color='b',width=.005,linewidth=.1,pivot='mid') # !!! pivot='mid','tip'

# deco
plt.plot()
plt.axis([0, Lx, 0, Ly])
plt.axis('equal')
plt.xlabel(r'$x$',fontsize=25)
plt.ylabel(r'$y$',fontsize=25)
theTitle = plt.title("Solution at time t=" +'{:04.2f}'.format(0*dt),
                     horizontalalignment='center',fontsize=20)
# circle
L=5
t=np.linspace(0,2*np.pi,300)
plt.plot(L+L*np.cos(t),L+L*np.sin(t),'r-',linewidth=3)

plt.draw()

#- B)  boucle
for iTime in range(0,nTime+1,jumpPrint):
    # load
    #extension = 'seed' + str(nSeed) + '_' + str(iTime).zfill(9) + '.udat'
    extension = str(iTime).zfill(9) + '.udat'
    x     = loadBinary('../data/particleX_'     + extension)
    y     = loadBinary('../data/particleY_'     + extension)
    theta = loadBinary('../data/particleTheta_' + extension)
    # update
    pltParticle.set_offsets(np.c_[x,y])
    pltParticle.set_UVC(lengthArrow*np.cos(theta),lengthArrow*np.sin(theta))
    # deco
    theTitle.set_text('Solution at time t=' +'{:04.2f}'.format(iTime*dt))
    plt.draw()
    if (saveVideo):
        plt.savefig('images/foo_'+str(iTime).zfill(9) +'.jpg', bbox_inches='tight')

sys.exit(0)




##---------------------------------------------------------------##
##---------------------------------------------------------------##

## mencoder ''mf://*.jpg'' -mf fps=40 -o ../bidon3.avi -ovc lavc -lavcopts vcodec=mpeg4
## mencoder ''mf://*.jpg'' -mf fps=25 -o ../bidon3_bis.avi -ovc x264 -lavcopts vcodec=mpeg4
## mencoder ''mf://*.jpg'' -mf fps=40 -o ../output3.avi -oac copy -ovc copy

