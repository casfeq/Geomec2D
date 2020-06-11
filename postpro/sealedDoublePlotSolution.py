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

# Get grid type for this run
gridType=[]
gridType.append("staggered")
# gridType.append("collocated+CDS")
# gridType.append("collocated+I2DPIS")

# Get dt
fileName=str(parentDirectory)+"/export/solveSealedDouble_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)
# dt=dt[0]

# Get timesteps exported
timesteps=[]
# timesteps.append("1")
# timesteps.append("31")
# timesteps.append("62")
timesteps.append("500")

# Define markers
markers=[]
markers.append("ko")
markers.append("ks")
markers.append("kx")

# Define labels
labels=[]
labels.append("Staggered")
# labels.append("Collocated (Non-stabilized)")
# labels.append("Collocated (Stabilized)")

# SUBPLOT

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,5))
fig.subplots_adjust(top=0.88,bottom=0.15,left=0.08,right=0.92,wspace=0.4)

# Define figure's name
plotName="plot/sealedDouble.png"

# Add subplot for pressure
fig.add_subplot(1,2,1)

# Plot pore-pressure
fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_PExact_"+gridType[0]+ \
	"-grid.txt"
pExact=np.loadtxt(fname=fileName)
pExact[:]=[x/1000 for x in pExact]
fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_YExact_"+gridType[0]+ \
	"-grid.txt"
yExact=np.loadtxt(fname=fileName)
exact,=plt.plot(pExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25,label="Analytical")
for j in range(0,len(gridType)):
	for i in range(0,len(timesteps)):
		fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_PPoreNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_YPNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		yNumeric=np.loadtxt(fname=fileName)
		if i==0:
			numeric,=plt.plot(pNumeric,yNumeric,markers[j],fillstyle='none',ms=5,mec='k', \
				mew=0.75,label=labels[j])
		else:
			numeric,=plt.plot(pNumeric,yNumeric,markers[j],fillstyle='none',ms=5,mec='k',mew=0.75)

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Add notes
# plt.text(0,3,"33.11 hours")
# plt.plot([2.6,7.5],[3.2,5],'-k',linewidth=0.5)
# plt.text(0,2.5,"85.54 days")
# plt.plot([2.6,6],[2.7,4],'-k',linewidth=0.5)
# plt.text(0,2,"344.93 days")
# plt.plot([2.6,5],[2.2,3],'-k',linewidth=0.5)
# plt.text(0,1.5,"689.85 days")
# plt.plot([2.6,4.2],[1.7,2],'-k',linewidth=0.5)

# Set axes' labels
plt.xlabel('Micropores Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add subplot for pressure
fig.add_subplot(1,2,2)

# Plot frac-pressure
fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_PExact_"+gridType[0]+ \
	"-grid.txt"
pExact=np.loadtxt(fname=fileName)
pExact[:]=[x/1000 for x in pExact]
fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_YExact_"+gridType[0]+ \
	"-grid.txt"
yExact=np.loadtxt(fname=fileName)
exact,=plt.plot(pExact,yExact,'-',color='k',fillstyle='none',linewidth=1.25)
for j in range(0,len(gridType)):
	for i in range(0,len(timesteps)):
		fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_PFracNumeric_dt=" \
			+format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/sealedDouble_"+solvedPairs[0]+"_YPNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		yNumeric=np.loadtxt(fname=fileName)
		numeric,=plt.plot(pNumeric,yNumeric,markers[j],fillstyle='none',ms=5,mec='k',mew=0.75)

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Add notes
# plt.text(0,3,"33.11 hours")
# plt.plot([2.6,7.5],[3.2,5],'-k',linewidth=0.5)
# plt.text(0,2.5,"85.54 days")
# plt.plot([2.6,6],[2.7,4],'-k',linewidth=0.5)
# plt.text(0,2,"344.93 days")
# plt.plot([2.6,5],[2.2,3],'-k',linewidth=0.5)
# plt.text(0,1.5,"689.85 days")
# plt.plot([2.6,4.2],[1.7,2],'-k',linewidth=0.5)

# Set axes' labels
plt.xlabel('Macropores Pressure (kPa)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add figure's legend
fig.legend(loc='upper center',ncol=2)

# Save figure
plt.savefig(plotName)

print("Plotted Sealed Column (Double Porosity)")