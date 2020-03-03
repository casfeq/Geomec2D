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

# Define grid type and interpolation scheme
gridType=[]
gridType.append("staggered")
gridType.append("collocated+CDS")
gridType.append("collocated+1DPIS")
gridType.append("collocated+I2DPIS")
gridType.append("collocated+C2DPIS")

# Define dt
dt=[]
dt.append(0.042726)
dt.append(0.017090)
dt.append(0.008545)
dt.append(0.001709)

# Define marker colors
colors=[]
colors.append("#fd411e")
colors.append("#f075e6")
colors.append("#0d75f8")
colors.append("#02c14d")
colors.append("#fec615")

# Define label
labels=[]
labels.append("Staggered grid")
labels.append("Collocated grid with CDS")
labels.append("Collocated grid with 1D-PIS")
labels.append("Collocated grid with I-2D-PIS")
labels.append("Collocated grid with C-2D-PIS")

# Define figure's name
plotName="plot/mandelStability_comparison.pdf"

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
		format(dt[i],".6f")+"_timeStep=1_staggered-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pExact[:]=[x/1000 for x in pExact]
	fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_staggered-grid.txt"
	xExact=np.loadtxt(fname=fileName)
	if i==0:
		exact,=plt.plot(xExact,pExact,'-k',fillstyle='none',label="Analytical")
	else:
		exact,=plt.plot(xExact,pExact,'-k',fillstyle='none')
	for j in range(0,len(gridType)):
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XPNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		xNumeric=np.loadtxt(fname=fileName)
		if i==0:
			numeric,=plt.plot(xNumeric,pNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5, \
				label=labels[j])
		numeric,=plt.plot(xNumeric,pNumeric,':.',color=colors[j],ms=10,mec='k',mew=0.5)

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

print("Plotted Mandel")