"""
	This source code is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Poroelasticity Problems". The routine here defined
	is used to plot the error obtained on the solution of the problem presented and solved by
	Terzaghi [1] against the meshsize to test the convergence of the method.

	Written by FERREIRA, C. A. S.

	Florianópolis, 2019.

	[1] TERZAGHI, K. Erdbaumechanik auf Bodenphysikalischer Grundlage. Franz Deuticke, Leipzig,
	1925.
"""

import numpy as np
import pathlib
import matplotlib.pyplot as plt
import math

"""    PLOT RESULTS
  ----------------------------------------------------------------"""

# Get parent directory
parentDirectory=pathlib.Path(__file__).resolve().parents[1]

# Get run info
runInfo=np.genfromtxt(parentDirectory/"export/convergenceTerzaghiTRunInfo.txt",dtype='str')

# Get grid type for this run
gridType=runInfo[0]

# Get interpolation scheme for this run if not staggered grid
if gridType!="staggered":
	interpScheme=runInfo[1]
	gridType=gridType+"+"+interpScheme

# Get h exported
h=[]
for i in range(0,len(runInfo)-2):
	h.append(runInfo[i+2])

# Define colors
colors=[]
colors.append("#fd411e")
colors.append("#f075e6")
colors.append("#0d75f8")
colors.append("#02c14d")

# SUBPLOT T-CONVERGENCE

# Define figure's name
plotName="plot/terzaghiTConvergence_"+gridType+"-grid.pdf"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,5))
fig.subplots_adjust(top=0.88,bottom=0.15,left=0.08,right=0.92,wspace=0.4)

# Add subplot for pressure
fig.add_subplot(1,2,1)

# Plot pressure
for i in range(0,len(h)):
	fileName=str(parentDirectory)+"/export/terzaghiErrorNorm_h="+h[i]+"_"+gridType+"-grid.txt"
	timeStepSize=np.loadtxt(fname=fileName,usecols=(0))
	pressureError=np.loadtxt(fname=fileName,usecols=(1))
	plt.loglog(timeStepSize,pressureError,'--o',color=colors[i],ms=5,mec='k',mew=0.5, \
		label="h="+str('%.4f'%float(h[i]))+" m")

# Plot triangle
e1=1.05*pressureError[0]
h1=0.85*timeStepSize[0]
h2=0.35*timeStepSize[0]
e2=e1*((h2/h1))
horCathetusX=[]
horCathetusY=[]
horCathetusX.append(h1)
horCathetusX.append(h2)
horCathetusY.append(e1)
horCathetusY.append(e1)
plt.loglog(horCathetusX,horCathetusY,'-k')
verCathetusX=[]
verCathetusY=[]
verCathetusX.append(h2)
verCathetusX.append(h2)
verCathetusY.append(e1)
verCathetusY.append(e2)
plt.loglog(verCathetusX,verCathetusY,'-k')
hypotenuseX=[]
hypotenuseY=[]
hypotenuseX.append(h1)
hypotenuseX.append(h2)
hypotenuseY.append(e1)
hypotenuseY.append(e2)
plt.loglog(hypotenuseX,hypotenuseY,'-k')

# Set axes' scale and limits
axes=plt.gca()
axes.set_xscale("log", nonposx='clip')
axes.set_yscale("log", nonposy='clip')
axes.text((h1+h2*2)/3,e1*1.05,"1")
axes.text(h2*0.85,(e1+e2*2)/3,"1")

# Set axes' labels
plt.xlabel('Timestep Size (s)')
plt.ylabel('Pressure Error (Pa)')
plt.grid(which='both',axis='both')

# Add subplot for vertical displacement
fig.add_subplot(1,2,2)

# Plot displacement
for i in range(0,len(h)):
	fileName=str(parentDirectory)+"/export/terzaghiErrorNorm_h="+h[i]+"_"+gridType+"-grid.txt"
	timeStepSize=np.loadtxt(fname=fileName,usecols=(0))
	displacementError=np.loadtxt(fname=fileName,usecols=(2))
	plt.loglog(timeStepSize,displacementError,'--o',color=colors[i],ms=5,mec='k',mew=0.5)

# Plot triangle
e1=1.05*displacementError[0]
h1=0.85*timeStepSize[0]
h2=0.35*timeStepSize[0]
e2=e1*((h2/h1))
horCathetusX=[]
horCathetusY=[]
horCathetusX.append(h1)
horCathetusX.append(h2)
horCathetusY.append(e1)
horCathetusY.append(e1)
plt.loglog(horCathetusX,horCathetusY,'-k')
verCathetusX=[]
verCathetusY=[]
verCathetusX.append(h2)
verCathetusX.append(h2)
verCathetusY.append(e1)
verCathetusY.append(e2)
plt.loglog(verCathetusX,verCathetusY,'-k')
hypotenuseX=[]
hypotenuseY=[]
hypotenuseX.append(h1)
hypotenuseX.append(h2)
hypotenuseY.append(e1)
hypotenuseY.append(e2)
plt.loglog(hypotenuseX,hypotenuseY,'-k')

# Set axes' scale and limits
axes=plt.gca()
axes.set_xscale("log", nonposx='clip')
axes.set_yscale("log", nonposy='clip')
axes.text((h1+h2*2)/3,e1*1.05,"1")
axes.text(h2*0.85,(e1+e2*2)/3,"1")

# Set axes' labels
plt.xlabel('Timestep Size (s)')
plt.ylabel('Vertical displacement error (m)')
plt.grid(which='both',axis='both')

# Add figure's legend
fig.legend(loc='upper center',ncol=4)

# Save figure
plt.savefig(plotName)

print("Plotted convergence results")