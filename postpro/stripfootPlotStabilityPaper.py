import numpy as np
import pathlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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
runInfo=np.genfromtxt(parentDirectory/"export/solveStripfootRunInfo.txt",dtype='str')

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

# Get dt
fileName=str(parentDirectory)+"/export/solveStripfoot_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)
dt=dt[0:4]

# Loope for grid type and interpolation schemes
for j in range(0,len(gridType)):

	# Create and define figure's size and margins
	fig=plt.figure(figsize=(8,7))
	fig.subplots_adjust(top=0.96,bottom=0.07,left=0.06,right=0.96,wspace=0.3,hspace=0.3)
		
	# Define figure's name
	plotName="plot/stripfootStability_"+gridType[j]+"-grid.pdf"

	# Loop for adding subplots
	ax=[]
	for i in range(0,len(dt)):

		# Add subplot
		ax.append(fig.add_subplot(2,2,i+1))

		# Plot pressure
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_PNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+format(dt[i], \
			".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		xPosition=np.loadtxt(fname=fileName)
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+format(dt[i], \
			".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		yPosition=np.loadtxt(fname=fileName)
		XPosition,YPosition=np.meshgrid(xPosition,yPosition)
		im=plt.contourf(XPosition,YPosition,pNumeric,500,cmap=cm.inferno)

		# This is the fix for the white lines between contour levels
		for c in im.collections:
			c.set_edgecolor("face")

		# Add colorbar
		cbar=fig.colorbar(im)
		cbar.ax.set_ylabel('Pressure (kPa)',rotation=90)

		# Set axes' labels
		plt.xlabel('Length (m)')
		plt.ylabel('Height (m)')		

	# Subplots' titles
	ax[0].title.set_text('25% of consolidation time')
	ax[1].title.set_text('10% of consolidation time')
	ax[2].title.set_text('5% of consolidation time')
	ax[3].title.set_text('1% of consolidation time')

	# Save figure
	plt.savefig(plotName)

# Create and define figure's size and margins
fig=plt.figure(figsize=(7,8))
fig.subplots_adjust(top=0.90,bottom=0.08,left=0.1,right=0.97,wspace=0.25,hspace=0.35)
	
# Define figure's name
plotName="plot/stripfootStability_along_line.pdf"

ax=[]
# Loop for grid type and interpolation schemes
for i in range(len(dt)-2,len(dt)):
	# Add subplot for micropores (x=strip)
	ax.append(fig.add_subplot(2,2,i-1))

	# Loope for grid type and interpolation schemes
	for j in range(0,len(gridType)):
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_PNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		xPosition=np.loadtxt(fname=fileName)
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		yPosition=np.loadtxt(fname=fileName)
		strip=int(len(xPosition)/5)
		if i == len(dt)-2:
			numeric,=plt.plot(pNumeric[:,strip],yPosition,markers[j],linewidth=0.75, \
				fillstyle='none',ms=5,mec='k',mew=0.75,label=labels[j])
		else:
			numeric,=plt.plot(pNumeric[:,strip],yPosition,markers[j],linewidth=0.75, \
				fillstyle='none',ms=5,mec='k',mew=0.75)

		# Set axes' scale and limits
		axes=plt.gca()
		axes.set_ylim([0,None])

		# Set axes' labels
		plt.xlabel('Pressure (kPa)')
		plt.ylabel('Height (m)')
		plt.grid(which='major',axis='both')

	# Add subplot for micropores (y=strip)
	ax.append(fig.add_subplot(2,2,i+1))

	# Loope for grid type and interpolation schemes
	for j in range(0,len(gridType)):
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_PNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		xPosition=np.loadtxt(fname=fileName)
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		yPosition=np.loadtxt(fname=fileName)
		strip=int(len(yPosition)/5*4)
		numeric,=plt.plot(xPosition,pNumeric[strip,:],markers[j],linewidth=0.75,fillstyle='none', \
			ms=5,mec='k',mew=0.75)

		# Set axes' scale and limits
		axes=plt.gca()
		axes.set_xlim([0,None])

		# Set axes' labels
		plt.xlabel('Length (m)')
		plt.ylabel('Pressure (kPa)')
		plt.grid(which='major',axis='both')

# Subplots' titles
ax[0].title.set_text('5% of consolidation time (x=1 m)')
ax[1].title.set_text('5% of consolidation time (y=4 m)')
ax[2].title.set_text('1% of consolidation time (x=1 m)')
ax[3].title.set_text('1% of consolidation time (y=4 m)')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

print("Plotted Stripfoot")