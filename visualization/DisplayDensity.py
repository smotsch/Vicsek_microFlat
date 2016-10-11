#!/usr/bin/python3

import time
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
import csv


# Parameters visualisation
shouldPlot1D = False
shouldPlotEnd = False
saveVideo = True

# ---------------- 0.1) Read parameters ----------------- #
# ------------------------------------------------------- #
f = open('../bin/PARAMETER_MicroVic_flat.txt')
lines = f.readlines()
f.close()
N         = int(lines[5-1])
Lx        = float(lines[12-1])
Ly        = float(lines[13-1])
dt        = float(lines[18-1])
Time      = float(lines[19-1])
nSeed     = int(lines[27-1])
dx        = float(lines[31-1])
dxy       = float(lines[32-1])
jumpPrint = int(lines[35-1])

# ------------------ 0.2) Initialisation ----------------- #
# -------------------------------------------------------- #
nTime = int(Time/dt+.5)
if (shouldPlot1D):
    # 1D
    nX = int(Lx/dx+.5)
    intX = np.arange(nX+1)*dx

    def loadBinary(nameFile):
        # Read the binary data in the file 'nameFile' in 'l'.
        # We also need to precise the number of rows (numberRow).
        f = open(nameFile, "rb")
        f.seek(4, os.SEEK_SET)
        l = np.fromfile(f, dtype=np.float64)
        f.close()
        return l

else:
    # 2D
    nX = int(Lx/dxy+.5)
    nY = int(Ly/dxy+.5)
    intX = np.arange(nX+1)*dxy
    intY = np.arange(nY+1)*dxy
    xx, yy = np.meshgrid(intX, intY, sparse=False)
    rhoMean = N/(nX*nY)/dxy/dxy
    zMax = 4*rhoMean

    def loadBinary2D(nameFile, numberRow):
        # Read the binary data in the file 'nameFile' in 'l'.
        # We also need to precise the number of rows (numberRow).
        f = open(nameFile, "rb")
        f.seek(4, os.SEEK_SET)
        x = np.fromfile(f, dtype=np.float64)
        M = x.reshape(numberRow, x.size/numberRow, order='F')  # Fortran style
        f.close()
        return M

if (shouldPlotEnd):
    jumpPrint = nTime

# init
# ----
if (saveVideo):
    plt.ioff()
else:
    plt.ion()
# plt.show()
plt.clf()


# --------------------------------------------------------------- #
# ------------------------  la boucle  -------------------------- #
# --------------------------------------------------------------- #

if (shouldPlot1D):
    # ------------------------  1D  -------------------------- #
    # -------------------------------------------------------- #
    # A) first image
    # load
    extension = 'seed' + str(nSeed) + '_' + str(0).zfill(9) + '.udat'
    rho1Dx = loadBinary('../data/dens1Dx_' + extension)
    v1Dx   = loadBinary('../data/v1Dx_'    + extension)
    # plot
    lines = plt.plot(intX, rho1Dx, intX, v1Dx, linewidth=4)
    # deco
    plt.axis([0, Lx, -1.2, 2.4])
    plt.xlabel(r'$x$', fontsize=25)
    theTitle = plt.title('Time t=' + '{:04.2f}'.format(0*dt),
                         horizontalalignment='center', fontsize=20)
    plt.legend([r'$\rho$', r'$v$'], fontsize=25)
    plt.legend(prop={'size': 25}, loc='upper left', fancybox=True, shadow=True)
    plt.draw()
    # B)  boucle
    for iTime in range(0, nTime+1, jumpPrint):
        # load
        extension = 'seed' + str(nSeed) + '_' + str(iTime).zfill(9) + '.udat'
        rho1Dx = loadBinary('../data/dens1Dx_' + extension);
        v1Dx   = loadBinary('../data/v1Dx_'    + extension)
        # update
        for line, y in zip(lines, [rho1Dx, v1Dx]):
            line.set_ydata(y)
        # deco
        theTitle.set_text('Time t=' + '{:04.2f}'.format(iTime*dt))
        plt.draw()
    if (saveVideo):
        plt.savefig('images/foo_' + str(iTime).zfill(9) + '.jpg',
                    bbox_inches='tight')
else:
    # ------------------------  2D  -------------------------- #
    # -------------------------------------------------------- #
    # A) first image
    # load
    extension = str(0).zfill(9) + '.udat'
    rho = loadBinary2D('../data/rho2D_' + extension, nX+1)
    u2D = loadBinary2D('../data/u2D_'   + extension, nX+1)
    v2D = loadBinary2D('../data/v2D_'   + extension, nX+1)
    # plot
    pltRho = plt.pcolormesh(xx, yy, rho.T, vmin=0, vmax=zMax,
                            shading='gouraud')  # 'gourand' for PIC method
    pltUV  = plt.quiver(xx, yy, u2D.T, v2D.T,
                        pivot='mid', scale=20,
                        width=.005, linewidth=.025, alpha=.5)
    # deco
    pltRho.set_cmap('hot_r')
    colorbar = plt.colorbar(pltRho)
    colorbar.set_ticks(np.arange(0, zMax+.01, zMax/5))
    plt.xlabel(r'$x$', fontsize=20)
    plt.ylabel(r'$y$', fontsize=20)
    theTitle = plt.title("Time  t=" + '{:04.2f}'.format(0*dt),
                         horizontalalignment='center')
    plt.draw()
    # B)  boucle
    for iTime in range(0, nTime+1, jumpPrint):
        # load
        extension = str(iTime).zfill(9) + '.udat'
        rho = loadBinary2D('../data/rho2D_' + extension, nX+1)
        u2D = loadBinary2D('../data/u2D_'   + extension, nX+1)
        v2D = loadBinary2D('../data/v2D_'   + extension, nX+1)
        # update
        pltRho.set_array(rho.T.ravel())
        pltUV.set_UVC(u2D.T, v2D.T)
        # deco
        theTitle.set_text('Time t=' + '{:04.2f}'.format(iTime*dt))
        plt.draw()
        if (saveVideo):
            plt.savefig('images/rhoUV_' + str(iTime).zfill(9) + '.jpg',
                        bbox_inches='tight')

sys.exit(0)

# --------------------------------------------------------------- #
# --------------------------------------------------------------- #

#  mencoder ''mf://*.jpg'' -mf fps=25 -o ../bidon3.avi -ovc lavc -lavcopts vcodec=mpeg4
#  mencoder ''mf://*.jpg'' -mf fps=25 -o ../bidon3_bis.avi -ovc x264 -lavcopts vcodec=mpeg4
#  mencoder ''mf://*.jpg'' -mf fps=25 -o ../bidon3_bis.avi -ovc x264 -lavcopts vcodec=mpeg4 -vf scale=640:-2
#  mencoder ''mf://*.jpg'' -mf fps=25 -o ../output3.avi -oac copy -ovc copy



# save csv
#---------
iTime = 2500                    # 0,750,1500,2500
extension = 'seed' + str(nSeed) + '_' + str(iTime).zfill(9) + '.udat'
rho1Dx = loadBinary('../data/dens1Dx_' + extension);
v1Dx   = loadBinary('../data/v1Dx_'    + extension)
with open('data_rho_v_time' + str(iTime) + '.csv', 'w') as fp:
    a = csv.writer(fp, delimiter=',')
    a.writerows( np.c_[intX,rho1Dx,v1Dx] )

iTime = 2500
rho = loadBinary2D('../data/rho2D_' + str(iTime).zfill(9) + '.udat',nY);
u2D = loadBinary2D('../data/u2D_'   + str(iTime).zfill(9) + '.udat',nY);
v2D = loadBinary2D('../data/v2D_'   + str(iTime).zfill(9) + '.udat',nY);
with open('data_rho2D_time' + str(iTime) + '.csv', 'w') as fp:
    a = csv.writer(fp, delimiter=',')
    a.writerows(np.c_[X.flatten(), Y.flatten(),
                      rho.flatten(), u2D.flatten(), v2D.flatten()] )


# ressources
# http://stackoverflow.com/questions/29009743/using-set-array-with-pyplot-pcolormesh-ruins-figure
