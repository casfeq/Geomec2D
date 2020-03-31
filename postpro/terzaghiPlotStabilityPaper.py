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

# Get dt
fileName=str(parentDirectory)+"/export/solveTerzaghi_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)
dt=dt[0:4]

# Define grid type and interpolation scheme
gridType=[]
gridType.append("staggered")
gridType.append("collocated+CDS")
gridType.append("collocated+I2DPIS")

# Define markers
markers=[]
markers.append(":ko")
markers.append(":ks")
markers.append(":kx")

# Define labels
labels=[]
labels.append("Staggered")
labels.append("Collocated (Non-stabilized)")
labels.append("Collocated (Stabilized)")

# Define figure's name
plotName="plot/terzaghiStability_comparison.pdf"

# Create and define figure's size and margins
fig=plt.figure(figsize=(7,9))
fig.subplots_adjust(top=0.90,bottom=0.06,left=0.06,right=0.98,wspace=0.2,hspace=0.25)

# Loop for adding pressure subplots
ax=[]
for i in range(0,len(dt)):

	# Add subplot
	ax.append(fig.add_subplot(2,2,i+1))

	# Plot pressure
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PExact_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_staggered-grid.txt"
	pExact=np.loadtxt(fname=fileName)
	pExact[:]=[x/1000 for x in pExact]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YExact_dt="+ \
		format(dt[i],".6f")+"_timeStep=1_staggered-grid.txt"
	yExact=np.loadtxt(fname=fileName)
	if i==0:
		exact,=plt.plot(pExact,yExact,'-',color='grey',fillstyle='none',label="Analytical")
	else:
		exact,=plt.plot(pExact,yExact,'-',color='grey',fillstyle='none')
	for j in range(0,len(gridType)):
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		yNumeric=np.loadtxt(fname=fileName)
		if i==0:
			numeric,=plt.plot(pNumeric,yNumeric,markers[j],linewidth=0.75,fillstyle='none',ms=5, \
				mec='k',mew=0.75,label=labels[j])
		else:
			numeric,=plt.plot(pNumeric,yNumeric,markers[j],linewidth=0.75,fillstyle='none',ms=5, \
				mec='k',mew=0.75)			

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
fig.legend(loc='upper center',ncol=2)

# Save figure
plt.savefig(plotName)

print("Plotted Terzaghi")