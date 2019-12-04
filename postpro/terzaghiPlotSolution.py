"""
	This source code is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Poroelasticity Problems". The routine here defined
	is used to plot the numerical and analytical results obtained on the solution of the problem
	presented and solved by Terzaghi [1].

	Written by FERREIRA, C. A. S.

	Florianópolis, 2019.

	[1] TERZAGHI, K. Erdbaumechanik auf Bodenphysikalischer Grundlage. Franz Deuticke, Leipzig,
	1925.
"""

import numpy as np
import pathlib
import matplotlib.pyplot as plt
import sys

"""    PLOT RESULTS
  ----------------------------------------------------------------"""

# Get pairs solved
solvedPairs=[]
for eachArg in sys.argv:
	solvedPairs.append(eachArg)
solvedPairs.pop(0)

# Get parent directory
parentDirectory=pathlib.Path(__file__).resolve().parents[1]

# Get run info
runInfo=np.genfromtxt(parentDirectory/"export/solveTerzaghiRunInfo.txt",dtype='str')

# Get grid type for this run
gridType=runInfo[0]

# Get interpolation scheme for this run if not staggered grid
if gridType!="staggered":
	interpScheme=runInfo[1]
	gridType=gridType+"+"+interpScheme

# Get dt
fileName=str(parentDirectory)+"/export/solveTerzaghi_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)

# Get timesteps exported
timesteps=[]
for i in range(0,len(runInfo)-2):
	timesteps.append(int(runInfo[i+2]))

# Define colors
colors=[]
colors.append("#fd411e")
colors.append("#f075e6")
colors.append("#0d75f8")
colors.append("#02c14d")

# SUBPLOT

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,5))
fig.subplots_adjust(top=0.88,bottom=0.15,left=0.08,right=0.92,wspace=0.4)

# Define figure's name
plotName="plot/terzaghiSolution_"+gridType+"-grid.png"

# Add subplot for pressure
fig.add_subplot(1,2,1)

# Plot pressure
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pExact[:]=[x/1000 for x in pExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	if i==0:
		exact,=plt.plot(pExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25, \
			label="Analytical")
	else:
		exact,=plt.plot(pExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(pNumeric,yNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5, \
		label="Timestep "+str(timesteps[i]))

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add subplot for pressure
fig.add_subplot(1,2,2)

# Plot displacement
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_VExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	vExact=np.loadtxt(fname=fileName)
	vExact[:]=[x*1000 for x in vExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	exact=plt.plot(vExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_VNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	vNumeric=np.loadtxt(fname=fileName)
	vNumeric[:]=[x*1000 for x in vNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YVNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric=plt.plot(vNumeric,yNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5)

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Vertical Displacement (mm)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

# SINGLE-PLOT

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,8))

# Define figure's name
plotName="plot/terzaghiSolution_p_"+gridType+"-grid.png"

# Plot pressure
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pExact[:]=[x/1000 for x in pExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	if i==0:
		exact,=plt.plot(pExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25, \
			label="Analytical")
	else:
		exact,=plt.plot(pExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(pNumeric,yNumeric,'.',color=colors[i],ms=10.0,mec='k',mew=0.5, \
		label="Timestep "+str(timesteps[i]))

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,8))

# Define figure's name
plotName="plot/terzaghiSolution_v_"+gridType+"-grid.png"

# Plot displacement
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_VExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	vExact=np.loadtxt(fname=fileName)
	vExact[:]=[x*1000 for x in vExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	if i==0:
		exact=plt.plot(vExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25, \
			label="Analytical")
	else:
		exact=plt.plot(vExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_VNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	vNumeric=np.loadtxt(fname=fileName)
	vNumeric[:]=[x*1000 for x in vNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YVNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric=plt.plot(vNumeric,yNumeric,'.',color=colors[i],ms=10.0,mec='k',mew=0.5, \
		label="Timestep "+str(timesteps[i]))

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Vertical Displacement (mm)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

print("Plotted Terzaghi")