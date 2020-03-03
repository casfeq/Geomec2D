"""
	This source code is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The routine here 
	defined is used to plot the numerical and analytical results obtained on the solution of the 
	problems presented and solved by Mandel [1] with parameters such that the stability of the 
	solution is tested.

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

# Define marker colors
colors=[]
colors.append("#fd411e")
colors.append("#f075e6")
colors.append("#0d75f8")
colors.append("#02c14d")

# Define figure's name
plotName="plot/mandelStability_p_"+gridType+"-grid.png"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,9))
fig.subplots_adjust(top=0.90,bottom=0.05,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# Loop for adding pressure subplots
ax=[]
for i in range(0,len(dt)):

	# Add subplot
	ax.append(fig.add_subplot(2,2,i+1))

	# Plot pressure
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PExact_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pExact[:]=[x/1000 for x in pExact]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	xExact=np.loadtxt(fname=fileName)
	if i==0:
		exact,=plt.plot(xExact,pExact,'-k',fillstyle='none',label="Analytical")
	else:
		exact,=plt.plot(xExact,pExact,'-k',fillstyle='none')
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PNumeric_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XPNumeric_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	xNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(xNumeric,pNumeric,':.',color=colors[i],ms=10,mec='k',mew=0.5, \
		label="$\Delta$t="+str(format(dt[i],'.2e'))+" s")

	# Set axes' scale and limits
	axes=plt.gca()
	axes.set_xlim([0,None])
	axes.set_ylim([0,None])

	# Set axes' labels
	plt.xlabel('Lenght (m)')
	plt.ylabel('Pressure (kPa)')
	plt.grid(which='major',axis='both')

# Subplots' titles
ax[0].title.set_text('25% of consolidation time')
ax[1].title.set_text('10% of consolidation time')
ax[2].title.set_text('5% of consolidation time')
ax[3].title.set_text('1% of consolidation time')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

# # Define figure's name
# plotName="plot/mandelStability_e_"+gridType+"-grid.png"

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,9))
# fig.subplots_adjust(top=0.90,bottom=0.05,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# # Loop for adding volumetric strain subplots
# for i in range(0,len(dt)):

# 	# Add subplot
# 	ax.append(fig.add_subplot(2,2,i+1))

# 	# Plot volumetric strain
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_EExact_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	eExact=np.loadtxt(fname=fileName)
# 	eExact[:]=[x*1e6 for x in eExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(xExact,eExact,'-k',fillstyle='none',label="Analytical")
# 	else:
# 		exact,=plt.plot(xExact,eExact,'-k',fillstyle='none')
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_ENumeric_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	eNumeric=np.loadtxt(fname=fileName)
# 	eNumeric[:]=[x*1e6 for x in eNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XENumeric_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,eNumeric,':.',color=colors[i],ms=10,mec='k',mew=0.5, \
# 		label="$\Delta$t="+str(format(dt[i],'.2e'))+" s")

# 	# Set axes' scale and limits
# 	axes=plt.gca()
# 	axes.set_xlim([0,None])

# 	# Set axes' labels
# 	plt.xlabel('Lenght (m)')
# 	plt.ylabel('Volumetric Strain (x$10^{-6}$)')
# 	plt.grid(which='major',axis='both')

# # Subplots' titles
# ax[0].title.set_text('25% of consolidation time')
# ax[1].title.set_text('10% of consolidation time')
# ax[2].title.set_text('5% of consolidation time')
# ax[3].title.set_text('1% of consolidation time')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

# # Define figure's name
# plotName="plot/mandelStability_u_"+gridType+"-grid.png"

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,9))
# fig.subplots_adjust(top=0.90,bottom=0.05,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# # Loop for adding horizontal displacement subplots
# for i in range(0,len(dt)):

# 	# Add subplot
# 	ax.append(fig.add_subplot(2,2,i+1))

# 	# Plot horizontal displacement
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UExact_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	uExact=np.loadtxt(fname=fileName)
# 	uExact[:]=[x*1e3 for x in uExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	xExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(xExact,uExact,'-k',fillstyle='none',label="Analytical")
# 	else:
# 		exact,=plt.plot(xExact,uExact,'-k',fillstyle='none')
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UNumeric_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	uNumeric=np.loadtxt(fname=fileName)
# 	uNumeric[:]=[x*1e3 for x in uNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XUNumeric_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	xNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(xNumeric,uNumeric,':.',color=colors[i],ms=10,mec='k',mew=0.5, \
# 		label="$\Delta$t="+str(format(dt[i],'.2e'))+" s")

# 	# Set axes' scale and limits
# 	axes=plt.gca()
# 	axes.set_xlim([0,None])

# 	# Set axes' labels
# 	plt.xlabel('Lenght (m)')
# 	plt.ylabel('Horizontal Displacement (mm)')
# 	plt.grid(which='major',axis='both')

# # Subplots' titles
# ax[0].title.set_text('25% of consolidation time')
# ax[1].title.set_text('10% of consolidation time')
# ax[2].title.set_text('5% of consolidation time')
# ax[3].title.set_text('1% of consolidation time')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

# # Define figure's name
# plotName="plot/mandelStability_v_"+gridType+"-grid.png"

# # Create and define figure's size and margins
# fig=plt.figure(figsize=(8,9))
# fig.subplots_adjust(top=0.90,bottom=0.05,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# # Loop for adding vertical displacement subplots
# for i in range(0,len(dt)):

# 	# Add subplot
# 	ax.append(fig.add_subplot(2,2,i+1))

# 	# Plot vertical displacement
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VExact_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	vExact=np.loadtxt(fname=fileName)
# 	vExact[:]=[x*1e3 for x in vExact]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YExact_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	yExact=np.loadtxt(fname=fileName)
# 	if i==0:
# 		exact,=plt.plot(vExact,yExact,'-k',fillstyle='none',label="Analytical")
# 	else:
# 		exact,=plt.plot(vExact,yExact,'-k',fillstyle='none')
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VNumeric_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	vNumeric=np.loadtxt(fname=fileName)
# 	vNumeric[:]=[x*1e3 for x in vNumeric]
# 	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YVNumeric_dt="+ \
# 		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
# 	yNumeric=np.loadtxt(fname=fileName)
# 	numeric,=plt.plot(vNumeric,yNumeric,':.',color=colors[i],ms=10,mec='k',mew=0.5, \
# 		label="$\Delta$t="+str(format(dt[i],'.2e'))+" s")

# 	# Set axes' scale and limits
# 	axes=plt.gca()
# 	axes.set_ylim([0,None])

# 	# Set axes' labels
# 	plt.xlabel('Vertical Displacement (mm)')
# 	plt.ylabel('Height (m)')
# 	plt.grid(which='major',axis='both')

# # Subplots' titles
# ax[0].title.set_text('25% of consolidation time')
# ax[1].title.set_text('10% of consolidation time')
# ax[2].title.set_text('5% of consolidation time')
# ax[3].title.set_text('1% of consolidation time')

# # Add figure's legend
# fig.legend(loc='upper center',ncol=3)

# # Save figure
# plt.savefig(plotName)

print("Plotted Mandel")