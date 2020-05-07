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

# Get dt
fileName=str(parentDirectory)+"/export/solveStripfoot_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)
dt=dt[0:2]

# Loope for grid type and interpolation schemes
for j in range(0,len(gridType)):

	# Create and define figure's size and margins
	fig=plt.figure(figsize=(8,7))
	fig.subplots_adjust(top=0.96,bottom=0.07,left=0.06,right=0.92,wspace=0.3,hspace=0.35)
		
	# Define figure's name
	plotName="plot/stripfootStability_"+gridType[j]+"-grid.png"

	# Loop for adding subplots
	ax=[]
	for i in range(0,len(dt)):

		# Add subplot for micropores
		ax.append(fig.add_subplot(2,2,2*i+1))

		# Plot pressure
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
		XPosition,YPosition=np.meshgrid(xPosition,yPosition)
		im=plt.contourf(XPosition,YPosition,pNumeric,500,cmap=cm.inferno)

		# Add colorbar
		cbar=fig.colorbar(im)
		cbar.ax.set_ylabel('Pressure (kPa)',rotation=90)

		# Set axes' labels
		plt.xlabel('Length (m)')
		plt.ylabel('Height (m)')

		# Add subplot for macropores
		ax.append(fig.add_subplot(2,2,2*i+2))

		# Plot pressure
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_MacroPNumeric_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		xPosition=np.loadtxt(fname=fileName)
		fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+ \
			format(dt[i],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
		yPosition=np.loadtxt(fname=fileName)
		XPosition,YPosition=np.meshgrid(xPosition,yPosition)
		im=plt.contourf(XPosition,YPosition,pNumeric,500,cmap=cm.inferno)

		# Add colorbar
		cbar=fig.colorbar(im)
		cbar.ax.set_ylabel('Pressure (kPa)',rotation=90)

		# Set axes' labels
		plt.xlabel('Length (m)')
		plt.ylabel('Height (m)')	

	# Subplots' titles
	ax[0].title.set_text('Micropores ($\Delta t=$'+str('%.2f'%float(dt[0]))+' s )')
	ax[1].title.set_text('Macropores ($\Delta t=$'+str('%.2f'%float(dt[0]))+' s )')
	ax[2].title.set_text('Micropores ($\Delta t=$'+str('%.2f'%float(dt[1]))+' s )')
	ax[3].title.set_text('Macropores ($\Delta t=$'+str('%.2f'%float(dt[1]))+' s )')

	# Save figure
	plt.savefig(plotName)

print("Plotted Stripfoot")