"""
	This source code is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The routine here 
	defined is used to plot the numerical and analytical results obtained on the solution of the 
	problem presented and solved by Mandel [1].

	Written by FERREIRA, C. A. S.

	Florianópolis, 2019.

	[1] MANDEL, J. Consolidation Des Sols (Étude Mathématique). Géotechnique, v. 3, n. 7, pp. 287-
	299, 1953.
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
runInfo=np.genfromtxt(parentDirectory/"export/solveMandelRunInfo.txt",dtype='str')

# Get grid type for this run
gridType=runInfo[0]

# Get interpolation scheme for this run if not staggered grid
if gridType!="staggered":
	interpScheme=runInfo[1]
	gridType=gridType+"+"+interpScheme

# Get dt
fileName=str(parentDirectory)+"/export/solveMandel_"+solvedPairs[0]+"RunInfo.txt"
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

# # SINGLE-PLOT

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,8))

# # Define figure's name
# plotName="plot/mandelSolution_p_"+gridType+"-grid.pdf"

# # Plot pressure
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	pExact=np.loadtxt(fname=fileName)
# 	pExact[:]=[x/1000 for x in pExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(xExact,pExact,'-',color='k',fillstyle='none',linewidth=1.25, \
# 			label="Analytical")
# 	else:
# 		exact,=plt.plot(xExact,pExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	pNumeric=np.loadtxt(fname=fileName)
# 	pNumeric[:]=[x/1000 for x in pNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XPNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,pNumeric,'.',color=colors[i],ms=10,mec='k',mew=0.5, \
# 		label="Time level "+str(timesteps[i]))

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_xlim([0,None])

# # Set axes' labels
# plt.xlabel('Length (m)')
# plt.ylabel('Pressure (kPa)')
# plt.grid(which='major',axis='both')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,8))

# # Define figure's name
# plotName="plot/mandelSolution_eps_"+gridType+"-grid.pdf"

# # Plot strain
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_EExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	epsExact=np.loadtxt(fname=fileName)
# 	epsExact[:]=[x*1000 for x in epsExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(xExact,epsExact,'-',color='k',fillstyle='none',linewidth=1.25, \
# 			label="Analytical")
# 	else:
# 		exact,=plt.plot(xExact,epsExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_ENumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	eNumeric=np.loadtxt(fname=fileName)
# 	eNumeric[:]=[x*1000 for x in eNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XENumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,eNumeric,'.',color=colors[i],ms=10,mec='k',mew=0.5, \
# 		label="Time level "+str(timesteps[i]))

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_xlim([0,None])

# # Set axes' labels
# plt.xlabel('Length (m)')
# plt.ylabel('Volumetric Strain (x$10^{-3}$)')
# plt.grid(which='major',axis='both')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,8))

# # Define figure's name
# plotName="plot/mandelSolution_u_"+gridType+"-grid.pdf"

# # Plot horizontal displacement
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	uExact=np.loadtxt(fname=fileName)
# 	uExact[:]=[x*1000 for x in uExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(xExact,uExact,'-',color='k',fillstyle='none',linewidth=1.25, \
# 			label="Analytical")
# 	else:
# 		exact,=plt.plot(xExact,uExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	uNumeric=np.loadtxt(fname=fileName)
# 	uNumeric[:]=[x*1000 for x in uNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XUNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,uNumeric,'.',color=colors[i],ms=10,mec='k',mew=0.5, \
# 		label="Time level "+str(timesteps[i]))

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_ylim([0,None])

# # Set axes' labels
# plt.xlabel('Length (m)')
# plt.ylabel('Horizontal Displacement (mm)')
# plt.grid(which='major',axis='both')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,8))

# # Define figure's name
# plotName="plot/mandelSolution_v_"+gridType+"-grid.pdf"

# # Plot vertical displacement
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	vExact=np.loadtxt(fname=fileName)
# 	vExact[:]=[x*1000 for x in vExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	yExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact=plt.plot(vExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25, \
# 			label="Analytical")
# 	else:
# 		exact=plt.plot(vExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	vNumeric=np.loadtxt(fname=fileName)
# 	vNumeric[:]=[x*1000 for x in vNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YVNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	yNumeric=np.loadtxt(fname=fileName)
# 	numeric=plt.plot(vNumeric,yNumeric,'.',color=colors[i],ms=10.0,mec='k',mew=0.5, \
# 		label="Time level "+str(timesteps[i]))

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_ylim([0,None])

# # Set axes' labels
# plt.xlabel('Vertical Displacement (mm)')
# plt.ylabel('Height (m)')
# plt.grid(which='major',axis='both')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

# SUBPLOT

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,8))
fig.subplots_adjust(top=0.88,bottom=0.15,left=0.08,right=0.92,wspace=0.4)

# Define figure's name
plotName="plot/mandelSolution_"+gridType+"-grid.pdf"

# Add subplot for pressure
fig.add_subplot(2,2,1)

# Plot pressure
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pExact[:]=[x/1000 for x in pExact]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	xExact=np.loadtxt(fname=fileName)
	if i==0:
		exact,=plt.plot(xExact,pExact,'-',color='k',fillstyle='none',linewidth=1.25, \
			label="Analytical")
	else:
		exact,=plt.plot(xExact,pExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XPNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	xNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(xNumeric,pNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5, \
		label="Time level "+str(timesteps[i]))

# Set axes' scale and limits
axes=plt.gca()
axes.set_xlim([0,None])

# Set axes' labels
plt.xlabel('Length (m)')
plt.ylabel('Pressure (kPa)')
plt.grid(which='major',axis='both')

# Add subplot for volumetric strain
fig.add_subplot(2,2,2)

# Plot strain
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_EExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	epsExact=np.loadtxt(fname=fileName)
	epsExact[:]=[x*1000 for x in epsExact]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	xExact=np.loadtxt(fname=fileName)
	exact,=plt.plot(xExact,epsExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_ENumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	eNumeric=np.loadtxt(fname=fileName)
	eNumeric[:]=[x*1000 for x in eNumeric]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XENumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	xNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(xNumeric,eNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5)

# Set axes' scale and limits
axes=plt.gca()
axes.set_xlim([0,None])

# Set axes' labels
plt.xlabel('Length (m)')
plt.ylabel('Volumetric Strain (x$10^{-3}$)')
plt.grid(which='major',axis='both')

# Add subplot for horizontal displacement
fig.add_subplot(2,2,3)

# Plot horizontal displacement
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	uExact=np.loadtxt(fname=fileName)
	uExact[:]=[x*1000 for x in uExact]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	xExact=np.loadtxt(fname=fileName)
	exact,=plt.plot(xExact,uExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	uNumeric=np.loadtxt(fname=fileName)
	uNumeric[:]=[x*1000 for x in uNumeric]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XUNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	xNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(xNumeric,uNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5)

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Length (m)')
plt.ylabel('Horizontal Displacement (mm)')
plt.grid(which='major',axis='both')

# Add subplot for vertical displacement
fig.add_subplot(2,2,4)

# Plot vertical displacement
for i in range(0,len(timesteps)):
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	vExact=np.loadtxt(fname=fileName)
	vExact[:]=[x*1000 for x in vExact]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YExact_dt="+format(dt,".6f")+ \
		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	exact=plt.plot(vExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VNumeric_dt="+ \
		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
	vNumeric=np.loadtxt(fname=fileName)
	vNumeric[:]=[x*1000 for x in vNumeric]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YVNumeric_dt="+ \
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

# # SUBPLOT 1

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,4))
# fig.subplots_adjust(top=0.85,bottom=0.15,left=0.08,right=0.92,wspace=0.4)

# # Define figure's name
# plotName="plot/mandelSolution_"+gridType+"-grid_1.pdf"

# # Add subplot for pressure
# fig.add_subplot(1,2,1)

# # Plot pressure
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	pExact=np.loadtxt(fname=fileName)
# 	pExact[:]=[x/1000 for x in pExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(xExact,pExact,'-',color='k',fillstyle='none',linewidth=1.25, \
# 			label="Analytical")
# 	else:
# 		exact,=plt.plot(xExact,pExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	pNumeric=np.loadtxt(fname=fileName)
# 	pNumeric[:]=[x/1000 for x in pNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XPNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,pNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5, \
# 		label="Time level "+str(timesteps[i]))

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_xlim([0,None])

# # Set axes' labels
# plt.xlabel('Length (m)')
# plt.ylabel('Pressure (kPa)')
# plt.grid(which='major',axis='both')

# # Add subplot for volumetric strain
# fig.add_subplot(1,2,2)

# # Plot strain
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_EExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	epsExact=np.loadtxt(fname=fileName)
# 	epsExact[:]=[x*1000 for x in epsExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	exact,=plt.plot(xExact,epsExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_ENumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	eNumeric=np.loadtxt(fname=fileName)
# 	eNumeric[:]=[x*1000 for x in eNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XENumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,eNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5)

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_xlim([0,None])

# # Set axes' labels
# plt.xlabel('Length (m)')
# plt.ylabel('Volumetric Strain (x$10^{-3}$)')
# plt.grid(which='major',axis='both')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

# # SUBPLOT 2

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,4))
# fig.subplots_adjust(top=0.85,bottom=0.15,left=0.08,right=0.92,wspace=0.4)

# # Define figure's name
# plotName="plot/mandelSolution_"+gridType+"-grid_2.pdf"

# # Add subplot for horizontal displacement
# fig.add_subplot(1,2,1)

# # Plot horizontal displacement
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	uExact=np.loadtxt(fname=fileName)
# 	uExact[:]=[x*1000 for x in uExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(xExact,uExact,'-',color='k',fillstyle='none',linewidth=1.25, \
# 			label="Analytical")
# 	else:
# 		exact,=plt.plot(xExact,uExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	uNumeric=np.loadtxt(fname=fileName)
# 	uNumeric[:]=[x*1000 for x in uNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XUNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,uNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5, \
# 		label="Time level "+str(timesteps[i]))

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_ylim([0,None])

# # Set axes' labels
# plt.xlabel('Length (m)')
# plt.ylabel('Horizontal Displacement (mm)')
# plt.grid(which='major',axis='both')

# # Add subplot for vertical displacement
# fig.add_subplot(1,2,2)

# # Plot vertical displacement
# for i in range(0,len(timesteps)):
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	vExact=np.loadtxt(fname=fileName)
# 	vExact[:]=[x*1000 for x in vExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YExact_dt="+format(dt,".6f")+ \
# 		"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	yExact=np.loadtxt(fname=fileName)
# 	exact=plt.plot(vExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	vNumeric=np.loadtxt(fname=fileName)
# 	vNumeric[:]=[x*1000 for x in vNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YVNumeric_dt="+ \
# 		format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType+"-grid.txt"
# 	yNumeric=np.loadtxt(fname=fileName)
# 	numeric=plt.plot(vNumeric,yNumeric,'.',color=colors[i],ms=7.5,mec='k',mew=0.5)

# # Set axes' scale and limits
# axes=plt.gca()
# axes.set_ylim([0,None])

# # Set axes' labels
# plt.xlabel('Vertical Displacement (mm)')
# plt.ylabel('Height (m)')
# plt.grid(which='major',axis='both')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

print("Plotted Mandel")