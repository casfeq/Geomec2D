"""
	This source code is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Poroelasticity Problems". The routine here defined
	is used to plot the numerical and analytical results obtained on the solution of the problems
	presented and solved by Terzaghi [1] with parameters such that the stability of the solution
	is tested.

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

# Get dtg[]
fileName=str(parentDirectory)+"/export/solveTerzaghi_"+solvedPairs[0]+"RunInfo.txt"
dtg=np.loadtxt(fname=fileName)

# Define marker colors
colors=[]
colors.append("#fd411e")
colors.append("#f075e6")
colors.append("#0d75f8")
colors.append("#02c14d")
colors.append("#fec615")

# Define figure's name
plotName="plot/terzaghiDimless_p_"+gridType+"-grid.pdf"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,9))
fig.subplots_adjust(top=0.90,bottom=0.05,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# Loop for adding pressure subplots
ax=[]
for i in range(0,4):

	# Add subplot
	ax.append(fig.add_subplot(2,2,i+1))

	# Plot pressure (analytical)
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PExact_dt="+ \
		format(dtg[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pMax=max(pExact)
	pExact[:]=[x/pMax for x in pExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dtg[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	yExact[:]=[x/6 for x in yExact]
	if i==0:
		exact,=plt.plot(pExact,yExact,'-k',fillstyle='none',label="Analytical")
	else:
		exact,=plt.plot(pExact,yExact,'-k',fillstyle='none')

	# Loop for plotting numerical solutions
	dtl=[]
	for j in range(0,len(solvedPairs)):
		fileName=str(parentDirectory)+"/export/solveTerzaghi_"+solvedPairs[j]+"RunInfo.txt"
		dtl=np.loadtxt(fname=fileName)
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_PExact_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		pExact=np.loadtxt(fname=fileName)
		pMax=max(pExact)
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_PNumeric_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/pMax for x in pNumeric]
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_YPNumeric_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		yNumeric=np.loadtxt(fname=fileName)
		yNumeric[:]=[x/6 for x in yNumeric]
		myFile=open(str(parentDirectory)+"/input/"+solvedPairs[j]+".txt","r")
		myLabel=myFile.readline().strip()
		myFile.close()
		if i==0:
			numeric,=plt.plot(pNumeric,yNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5, \
				label=str(myLabel))
		else:
			numeric,=plt.plot(pNumeric,yNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5)

	# Set axes' scale and limits
	axes=plt.gca()
	axes.set_xlim([0,None])
	axes.set_ylim([0,None])

	# Set axes' labels
	plt.xlabel('Dimensionless Pressure')
	plt.ylabel('Dimensionless Height')
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

# Define figure's name
plotName="plot/terzaghiDimless_p_"+gridType+"-grid_1.pdf"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,4.5))
fig.subplots_adjust(top=0.80,bottom=0.10,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# Loop for adding pressure subplots
ax=[]
for i in range(0,2):

	# Add subplot
	ax.append(fig.add_subplot(1,2,i+1))

	# Plot pressure (analytical)
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PExact_dt="+ \
		format(dtg[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pMax=max(pExact)
	pExact[:]=[x/pMax for x in pExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dtg[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	yExact[:]=[x/6 for x in yExact]
	if i==0:
		exact,=plt.plot(pExact,yExact,'-k',fillstyle='none',label="Analytical")
	else:
		exact,=plt.plot(pExact,yExact,'-k',fillstyle='none')

	# Loop for plotting numerical solutions
	dtl=[]
	for j in range(0,len(solvedPairs)):
		fileName=str(parentDirectory)+"/export/solveTerzaghi_"+solvedPairs[j]+"RunInfo.txt"
		dtl=np.loadtxt(fname=fileName)
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_PExact_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		pExact=np.loadtxt(fname=fileName)
		pMax=max(pExact)
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_PNumeric_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/pMax for x in pNumeric]
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_YPNumeric_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		yNumeric=np.loadtxt(fname=fileName)
		yNumeric[:]=[x/6 for x in yNumeric]
		myFile=open(str(parentDirectory)+"/input/"+solvedPairs[j]+".txt","r")
		myLabel=myFile.readline().strip()
		myFile.close()
		if i==0:
			numeric,=plt.plot(pNumeric,yNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5, \
				label=str(myLabel))
		else:
			numeric,=plt.plot(pNumeric,yNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5)

	# Set axes' scale and limits
	axes=plt.gca()
	axes.set_xlim([0,None])
	axes.set_ylim([0,None])

	# Set axes' labels
	plt.xlabel('Dimensionless Pressure')
	plt.ylabel('Dimensionless Height')
	plt.grid(which='major',axis='both')

# Subplots' titles
ax[0].title.set_text('25% of consolidation time')
ax[1].title.set_text('10% of consolidation time')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

# Define figure's name
plotName="plot/terzaghiDimless_p_"+gridType+"-grid_2.pdf"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,4.5))
fig.subplots_adjust(top=0.80,bottom=0.10,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# Loop for adding pressure subplots
ax=[]
for i in range(2,4):

	# Add subplot
	ax.append(fig.add_subplot(1,2,i-1))

	# Plot pressure (analytical)
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PExact_dt="+ \
		format(dtg[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pMax=max(pExact)
	pExact[:]=[x/pMax for x in pExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dtg[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	yExact[:]=[x/6 for x in yExact]
	if i==2:
		exact,=plt.plot(pExact,yExact,'-k',fillstyle='none',label="Analytical")
	else:
		exact,=plt.plot(pExact,yExact,'-k',fillstyle='none')

	# Loop for plotting numerical solutions
	dtl=[]
	for j in range(0,len(solvedPairs)):
		fileName=str(parentDirectory)+"/export/solveTerzaghi_"+solvedPairs[j]+"RunInfo.txt"
		dtl=np.loadtxt(fname=fileName)
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_PExact_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		pExact=np.loadtxt(fname=fileName)
		pMax=max(pExact)
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_PNumeric_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/pMax for x in pNumeric]
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[j]+"_YPNumeric_dt="+ \
			format(dtl[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
		yNumeric=np.loadtxt(fname=fileName)
		yNumeric[:]=[x/6 for x in yNumeric]
		myFile=open(str(parentDirectory)+"/input/"+solvedPairs[j]+".txt","r")
		myLabel=myFile.readline().strip()
		myFile.close()
		if i==2:
			numeric,=plt.plot(pNumeric,yNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5, \
				label=str(myLabel))
		else:
			numeric,=plt.plot(pNumeric,yNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5)

	# Set axes' scale and limits
	axes=plt.gca()
	axes.set_xlim([0,None])
	axes.set_ylim([0,None])

	# Set axes' labels
	plt.xlabel('Dimensionless Pressure')
	plt.ylabel('Dimensionless Height')
	plt.grid(which='major',axis='both')

# Subplots' titles
ax[0].title.set_text('5% of consolidation time')
ax[1].title.set_text('1% of consolidation time')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

print("Plotted Terzaghi")