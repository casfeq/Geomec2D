"""
	The routine here defined is used to plot the numerical and analytical results obtained on the
	solution of consolidation problems with parameters such that the stability of the solution is
	tested.

	Written by FERREIRA, C. A. S.

	Florianópolis, 2020.
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

# Define marker colors
colors=[]
colors.append("#fd411e")
colors.append("#f075e6")
colors.append("#0d75f8")
colors.append("#02c14d")

# Define figure's name
plotName="plot/terzaghiStability_p_"+gridType+"-grid.png"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,9))
fig.subplots_adjust(top=0.90,bottom=0.05,left=0.08,right=0.96,wspace=0.4,hspace=0.3)

# Loop for adding pressure subplots
ax=[]
for i in range(0,len(dt)):

	# Add subplot
	ax.append(fig.add_subplot(2,2,i+1))

	# Plot pressure
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PNumeric_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_"+gridType+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(pNumeric,yNumeric,':.',color=colors[i],ms=10,mec='k',mew=0.5, \
		label="$\Delta$t="+str(format(dt[i],'.2e'))+" s")

	# Set axes' scale and limits
	axes=plt.gca()
	axes.set_xlim([0,None])
	axes.set_ylim([0,None])

	# Set axes' labels
	plt.xlabel('Pressure (kPa)')
	plt.ylabel('Height (m)')
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

print("Plotted Terzaghi")