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

# Get dt
fileName=str(parentDirectory)+"/export/solveTerzaghi_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)

# Define plot style
plotstyle=[]
plotstyle.append(":s")
plotstyle.append(":o")
plotstyle.append(":^")
plotstyle.append(":d")
plotstyle.append(":*")
plotstyle.append(":x")

# Define marker colors
colors=[]
colors.append("#414042") # dark gray
colors.append("#fd411e") # red
# colors.append("#f075e6") # magenta
colors.append("#0d75f8") # blue
colors.append("#02c14d") # green

# Define grids used
gridType=[]
gridType.append("staggered")
gridType.append("collocated+CDS")
gridType.append("collocated+I2DPIS")

# Define figure's name
plotName="plot/doublePorosity_dt_stability.png"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,9))
fig.subplots_adjust(top=0.93,bottom=0.05,left=0.08,right=0.96,wspace=0.2,hspace=0.25)

# Add subplot 1 for micro-pressure
fig.add_subplot(2,2,2)
plt.subplot(2,2,2).set_title("Micro-pores ($\Delta t=$"+str(format(dt[0],'.2e'))+" s)",fontsize=10)

# Plot micro-pressure 1
for j in range(0,len(gridType)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PNumeric_dt="+format(dt[0], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
		format(dt[0],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(pNumeric,yNumeric,plotstyle[j],color=colors[j],ms=6,mec='k',mew=0.25, \
		label=gridType[j])

# Set axes' scale and limits
axes=plt.gca()
# axes.set_xlim([0,None])
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add subplot 2 for micro-pressure
fig.add_subplot(2,2,4)
plt.subplot(2,2,4).set_title("Micro-pores ($\Delta t=$"+str(format(dt[1],'.2e'))+" s)",fontsize=10)

# Plot micro-pressure 2
for j in range(0,len(gridType)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_PNumeric_dt="+format(dt[1], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
		format(dt[1],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(pNumeric,yNumeric,plotstyle[j],color=colors[j],ms=6,mec='k',mew=0.25)

# Set axes' scale and limits
axes=plt.gca()
# axes.set_xlim([0,None])
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add subplot 1 for macro-pressure
fig.add_subplot(2,2,1)
plt.subplot(2,2,1).set_title("Macro-pores ($\Delta t=$"+str(format(dt[0],'.2e'))+" s)",fontsize=10)

# Plot macro-pressure 1
for j in range(0,len(gridType)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_MacroPNumeric_dt="+ \
		format(dt[0],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
		format(dt[0],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(pNumeric,yNumeric,plotstyle[j],color=colors[j],ms=6,mec='k',mew=0.25)

# Set axes' scale and limits
axes=plt.gca()
axes.set_xlim([0,None])
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add subplot 2 for macro-pressure
fig.add_subplot(2,2,3)
plt.subplot(2,2,3).set_title("Macro-pores ($\Delta t=$"+str(format(dt[1],'.2e'))+" s)",fontsize=10)

# Plot macro-pressure 2
for j in range(0,len(gridType)):
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_MacroPNumeric_dt="+ \
		format(dt[1],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/terzaghi_"+solvedPairs[0]+"_YPNumeric_dt="+ \
		format(dt[1],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yNumeric=np.loadtxt(fname=fileName)
	numeric,=plt.plot(pNumeric,yNumeric,plotstyle[j],color=colors[j],ms=6,mec='k',mew=0.25)

# Set axes' scale and limits
axes=plt.gca()
axes.set_xlim([0,None])
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add figure's legend
fig.legend(loc='upper center',ncol=3)

# Save figure
plt.savefig(plotName)

print("Plotted Terzaghi")

# Get dt
fileName=str(parentDirectory)+"/export/solveStripfoot_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)

# Loop for stripfoot
for j in range(0,len(gridType)):

	# Define figure's name
	plotName="plot/doublePorosity_stripfoot_dt_stability_"+gridType[j]+"-grid.png"

	# Create and define figure's size and margins
	fig=plt.figure(figsize=(8,8))
	fig.subplots_adjust(top=0.97,bottom=0.05,left=0.08,right=0.96,wspace=0.2,hspace=0.25)

	# Add subplot 1 for micro-pressure
	fig.add_subplot(2,2,2)
	plt.subplot(2,2,2).set_title("Micro-pores ($\Delta t=$"+str(format(dt[0],'.2e'))+" s)", \
		fontsize=10)

	# Plot micro-pressure 1
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_PNumeric_dt="+ \
		format(dt[0],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+format(dt[0], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	xPosition=np.loadtxt(fname=fileName)
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+format(dt[0], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yPosition=np.loadtxt(fname=fileName)
	XPosition,YPosition=np.meshgrid(xPosition,yPosition)
	plt.contourf(XPosition,YPosition,pNumeric,100,cmap=cm.inferno)

	# Set axes' labels
	plt.xlabel('Length (m)')
	plt.ylabel('Height (m)')

	# Add subplot 2 for micro-pressure
	fig.add_subplot(2,2,4)
	plt.subplot(2,2,4).set_title("Micro-pores ($\Delta t=$"+str(format(dt[1],'.2e'))+" s)", \
		fontsize=10)

	# Plot micro-pressure 2
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_PNumeric_dt="+ \
		format(dt[1],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+format(dt[1], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	xPosition=np.loadtxt(fname=fileName)
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+format(dt[1], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yPosition=np.loadtxt(fname=fileName)
	XPosition,YPosition=np.meshgrid(xPosition,yPosition)
	plt.contourf(XPosition,YPosition,pNumeric,100,cmap=cm.inferno)

	# Set axes' labels
	plt.xlabel('Length (m)')
	plt.ylabel('Height (m)')

	# Add subplot 1 for macro-pressure
	fig.add_subplot(2,2,1)
	plt.subplot(2,2,1).set_title("Micro-pores ($\Delta t=$"+str(format(dt[0],'.2e'))+" s)", \
		fontsize=10)

	# Plot macro-pressure 1
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_MacroPNumeric_dt="+ \
		format(dt[0],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+format(dt[0], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	xPosition=np.loadtxt(fname=fileName)
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+format(dt[0], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yPosition=np.loadtxt(fname=fileName)
	XPosition,YPosition=np.meshgrid(xPosition,yPosition)
	plt.contourf(XPosition,YPosition,pNumeric,100,cmap=cm.inferno)

	# Set axes' labels
	plt.xlabel('Length (m)')
	plt.ylabel('Height (m)')

	# Add subplot 2 for macro-pressure
	fig.add_subplot(2,2,3)
	plt.subplot(2,2,3).set_title("Micro-pores ($\Delta t=$"+str(format(dt[1],'.2e'))+" s)", \
		fontsize=10)

	# Plot macro-pressure 2
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_MacroPNumeric_dt="+ \
		format(dt[1],".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	pNumeric=np.loadtxt(fname=fileName)
	pNumeric[:]=[x/1000 for x in pNumeric]
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_xCoord_dt="+format(dt[1], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	xPosition=np.loadtxt(fname=fileName)
	fileName=str(parentDirectory)+"/export/stripfoot_"+solvedPairs[0]+"_yCoord_dt="+format(dt[1], \
		".6f")+"_timeStep=1_"+gridType[j]+"-grid.txt"
	yPosition=np.loadtxt(fname=fileName)
	XPosition,YPosition=np.meshgrid(xPosition,yPosition)
	plt.contourf(XPosition,YPosition,pNumeric,100,cmap=cm.inferno)

	# Set axes' labels
	plt.xlabel('Length (m)')
	plt.ylabel('Height (m)')

	# Save figure
	plt.savefig(plotName)

print("Plotted Stripfoot")