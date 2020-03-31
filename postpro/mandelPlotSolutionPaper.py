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
gridType=[]
gridType.append("staggered")
gridType.append("collocated+CDS")
gridType.append("collocated+I2DPIS")

# Get dt
fileName=str(parentDirectory)+"/export/solveMandel_"+solvedPairs[0]+"RunInfo.txt"
dt=np.loadtxt(fname=fileName)
dt=dt[0]

# Get timesteps exported
timesteps=[]
timesteps.append("1")
timesteps.append("31")
timesteps.append("125")
timesteps.append("500")

# Define markers
markers=[]
markers.append("ko")
markers.append("ks")
markers.append("kx")

# Define labels
labels=[]
labels.append("Staggered")
labels.append("Collocated (Non-stabilized)")
labels.append("Collocated (Stabilized)")

# SUBPLOT

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,7))
fig.subplots_adjust(top=0.90,bottom=0.065,left=0.08,right=0.98,wspace=0.4)

# Define figure's name
plotName="plot/mandelSolution_paper.pdf"

# Add subplot for pressure
fig.add_subplot(2,2,1)

# Plot pressure
for j in range(0,len(gridType)):
	for i in range(0,len(timesteps)):
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		pExact=np.loadtxt(fname=fileName)
		pExact[:]=[x/1000 for x in pExact]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		xExact=np.loadtxt(fname=fileName)
		if i==0 and j==0:
			exact,=plt.plot(xExact,pExact,'-',color='grey',fillstyle='none',linewidth=1.25, \
				label="Analytical")
		if j==0:
			exact,=plt.plot(xExact,pExact,'-',color='grey',fillstyle='none',linewidth=1.25)
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_PNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		pNumeric=np.loadtxt(fname=fileName)
		pNumeric[:]=[x/1000 for x in pNumeric]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XPNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		xNumeric=np.loadtxt(fname=fileName)
		if i==0:
			numeric,=plt.plot(xNumeric,pNumeric,markers[j],fillstyle='none',ms=5,mec='k',mew=0.75, \
				label=labels[j])
		else:
			numeric,=plt.plot(xNumeric,pNumeric,markers[j],fillstyle='none',ms=5,mec='k',mew=0.75)

# Set axes' scale and limits
axes=plt.gca()
axes.set_xlim([0,None])

# Set axes' labels
plt.xlabel('Length (m)')
plt.ylabel('Pressure (kPa)')
plt.grid(which='major',axis='both')

# Add notes
plt.text(0.05,20,"33.11 hours")
plt.plot([1.5,4.5],[22,40],'-k',linewidth=0.5)
plt.text(0.05,15,"42.77 days")
plt.plot([1.4,4],[17,28],'-k',linewidth=0.5)
plt.text(0.05,10,"172.46 days")
plt.plot([1.6,3.7],[12,20],'-k',linewidth=0.5)
plt.text(0.05,5,"689.85 days")
plt.plot([1.6,3.7],[7,10],'-k',linewidth=0.5)

# Add subplot for volumetric strain
fig.add_subplot(2,2,2)

# Plot strain
for j in range(0,len(gridType)):
	for i in range(0,len(timesteps)):
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_EExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		epsExact=np.loadtxt(fname=fileName)
		epsExact[:]=[x*1000 for x in epsExact]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		xExact=np.loadtxt(fname=fileName)
		if j==0:
			exact,=plt.plot(xExact,epsExact,'-',color='grey',fillstyle='none',linewidth=1.25)
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_ENumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		eNumeric=np.loadtxt(fname=fileName)
		eNumeric[:]=[x*1000 for x in eNumeric]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XENumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		xNumeric=np.loadtxt(fname=fileName)
		numeric,=plt.plot(xNumeric,eNumeric,markers[j],fillstyle='none',ms=5,mec='k',mew=0.75)

# Set axes' scale and limits
axes=plt.gca()
axes.set_xlim([0,None])

# Set axes' labels
plt.xlabel('Length (m)')
plt.ylabel('Volumetric Strain (x$10^{-3}$)')
plt.grid(which='major',axis='both')

# Add notes
plt.text(0.05,-0.0250,"33.11 hours")
plt.plot([1.5,4.5],[-0.0240,-0.0075],'-k',linewidth=0.5)
plt.text(0.05,-0.0275,"42.77 days")
plt.plot([1.4,4],[-0.0265,-0.0150],'-k',linewidth=0.5)
plt.text(0.05,-0.0300,"172.46 days")
plt.plot([1.6,4.5],[-0.0290,-0.0250],'-k',linewidth=0.5)
plt.text(0.05,-0.0325,"689.85 days")
plt.plot([1.6,4.5],[-0.0315,-0.0300],'-k',linewidth=0.5)

# Add subplot for horizontal displacement
fig.add_subplot(2,2,3)

# Plot horizontal displacement
for j in range(0,len(gridType)):
	for i in range(0,len(timesteps)):
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		uExact=np.loadtxt(fname=fileName)
		uExact[:]=[x*1000 for x in uExact]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		xExact=np.loadtxt(fname=fileName)
		if j==0:
			exact,=plt.plot(xExact,uExact,'-',color='grey',fillstyle='none',linewidth=1.25)
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_UNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		uNumeric=np.loadtxt(fname=fileName)
		uNumeric[:]=[x*1000 for x in uNumeric]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_XUNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		xNumeric=np.loadtxt(fname=fileName)
		numeric,=plt.plot(xNumeric,uNumeric,markers[j],fillstyle='none',ms=5,mec='k',mew=0.75)

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Length (m)')
plt.ylabel('Horizontal Displacement (mm)')
plt.grid(which='major',axis='both')

# Add notes
plt.text(0,0.14,"33.11 hours")
plt.plot([1.6,5],[0.142,0.147],'-k',linewidth=0.5)
plt.text(0,0.125,"42.77 days")
plt.plot([1.5,5],[0.128,0.139],'-k',linewidth=0.5)
plt.text(0,0.11,"172.46 days")
plt.plot([1.6,5],[0.115,0.126],'-k',linewidth=0.5)
plt.text(0,0.095,"689.85 days")
plt.plot([1.65,5],[0.0975,0.105],'-k',linewidth=0.5)

# Add subplot for vertical displacement
fig.add_subplot(2,2,4)

# Plot vertical displacement
for j in range(0,len(gridType)):
	for i in range(0,len(timesteps)):
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		vExact=np.loadtxt(fname=fileName)
		vExact[:]=[x*1000 for x in vExact]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YExact_dt="+format(dt,".6f")+ \
			"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		yExact=np.loadtxt(fname=fileName)
		if j==0:
			exact=plt.plot(vExact,yExact,'-',color='grey',fillstyle='none',linewidth=1.25)
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_VNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		vNumeric=np.loadtxt(fname=fileName)
		vNumeric[:]=[x*1000 for x in vNumeric]
		fileName=str(parentDirectory)+"/export/mandel_"+solvedPairs[0]+"_YVNumeric_dt="+ \
			format(dt,".6f")+"_timeStep="+str(timesteps[i])+"_"+gridType[j]+"-grid.txt"
		yNumeric=np.loadtxt(fname=fileName)
		numeric,=plt.plot(vNumeric,yNumeric,markers[j],fillstyle='none',ms=5,mec='k',mew=0.75)

# Set axes' scale and limits
axes=plt.gca()
axes.set_ylim([0,None])

# Set axes' labels
plt.xlabel('Vertical Displacement (mm)')
plt.ylabel('Height (m)')
plt.grid(which='major',axis='both')

# Add notes
plt.text(-0.075,4.5,"33.11 hours")
plt.plot([-0.165,-0.077],[4.6,4.6],'-k',linewidth=0.5)
plt.text(-0.075,4,"42.77 days")
plt.plot([-0.155,-0.077],[4.1,4.1],'-k',linewidth=0.5)
plt.text(-0.225,2,"172.46 days")
plt.plot([-0.175,-0.175],[2.3,4.4],'-k',linewidth=0.5)
plt.text(-0.225,1.5,"689.85 days")
plt.plot([-0.15,-0.10],[1.6,2.2],'-k',linewidth=0.5)

# Add figure's legend
fig.legend(loc='upper center',ncol=2)

# Save figure
plt.savefig(plotName)

print("Plotted Mandel")