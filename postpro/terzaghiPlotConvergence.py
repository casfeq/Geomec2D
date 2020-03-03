"""
	This source code is part of the development of a master's thesis entitled "Analysis of Numerical
	Schemes in Collocated and Staggered Grids for Problems of Poroelasticity". The routine here 
	defined is used to plot the error obtained on the solution of the problem presented and solved
	by Terzaghi [1] against the meshsize to test the convergence of the method.

	Written by FERREIRA, C. A. S.

	Florianópolis, 2019.

	[1] TERZAGHI, K. Erdbaumechanik auf Bodenphysikalischer Grundlage. Franz Deuticke, Leipzig,
	1925.
"""

import numpy as np
import pathlib
import matplotlib.pyplot as plt
import math

"""    PLOT RESULTS
  ----------------------------------------------------------------"""

# Get parent directory
parentDirectory=pathlib.Path(__file__).resolve().parents[1]

# Get run info
runInfo=np.genfromtxt(parentDirectory/"export/convergenceTerzaghiRunInfo.txt",dtype='str')

# Get grid type for this run
gridType=runInfo[0]

# Get interpolation scheme for this run if not staggered grid
if gridType!="staggered":
	interpScheme=runInfo[1]
	gridType=gridType+"+"+interpScheme

# Get dt exported
dt=[]
for i in range(0,len(runInfo)-2):
	dt.append(runInfo[i+2])

# Define colors
colors=[]
colors.append("#fd411e")
colors.append("#f075e6")
colors.append("#0d75f8")
colors.append("#02c14d")

# # Defines linear regression
# pointsUsed=3
# Coef=[]
# def lin(h):
# 	return 10**(Coef[0]+Coef[1]*math.log10(h))

# SUBPLOT H-CONVERGENCE

# Define figure's name
plotName="plot/terzaghiConvergence_"+gridType+"-grid.pdf"

# Create and define figure's size and margins
fig=plt.figure(figsize=(8,5))
fig.subplots_adjust(top=0.88,bottom=0.15,left=0.08,right=0.92,wspace=0.4)

# Add subplot for pressure
fig.add_subplot(1,2,1)

# Plot pressure
# xData=[]
# yData=[]
for i in range(0,len(dt)):
	fileName=str(parentDirectory)+"/export/terzaghiErrorNorm_dt="+dt[i]+"_"+gridType+"-grid.txt"
	charLength=np.loadtxt(fname=fileName,usecols=(0))
	pressureError=np.loadtxt(fname=fileName,usecols=(1))
	plt.loglog(charLength,pressureError,'--o',color=colors[i],ms=5,mec='k',mew=0.5, \
		label="$\Delta$t="+str('%.0f'%float(dt[i]))+" s")
	# data=np.loadtxt(fname=fileName)
	# for i in range(0,pointsUsed):
	# 	row=[]
	# 	row.append(1)
	# 	row.append(math.log10(data[i,0]))
	# 	xData.append(row)
	# 	yData.append(math.log10(data[i,1]))

# # Determine coefficients for linear regression
# xDataT=xData
# xData=np.transpose(xData)
# yData=np.transpose(yData)
# Coef=np.matmul(xData,xDataT)
# Coef=np.linalg.inv(Coef)
# yData=np.matmul(xData,yData)
# Coef=np.matmul(Coef,yData)

# # Plot linear regression
# hReg=[]
# linReg=[]
# for i in range(0,len(charLength)):
# 	hReg.append(charLength[i])
# 	linReg.append(lin(hReg[i]))
# plt.loglog(hReg,linReg,'--k',linewidth=1.2,label='Lin.Reg.('+str(round(Coef[1],1))+')')

# Plot triangle
e1=1.05*pressureError[0]
h1=0.95*charLength[0]
h2=0.95*charLength[2]
e2=e1*((h2/h1)**2)
horCathetusX=[]
horCathetusY=[]
horCathetusX.append(h1)
horCathetusX.append(h2)
horCathetusY.append(e1)
horCathetusY.append(e1)
plt.loglog(horCathetusX,horCathetusY,'-k')
verCathetusX=[]
verCathetusY=[]
verCathetusX.append(h2)
verCathetusX.append(h2)
verCathetusY.append(e1)
verCathetusY.append(e2)
plt.loglog(verCathetusX,verCathetusY,'-k')
hypotenuseX=[]
hypotenuseY=[]
hypotenuseX.append(h1)
hypotenuseX.append(h2)
hypotenuseY.append(e1)
hypotenuseY.append(e2)
plt.loglog(hypotenuseX,hypotenuseY,'-k')

# Set axes' scale and limits
axes=plt.gca()
axes.set_xscale("log", nonposx='clip')
axes.set_yscale("log", nonposy='clip')
axes.text((h1+h2*2)/3,e1*1.05,"1")
axes.text(h2*0.9,(e1+e2*2)/3,"2")

# Set axes' labels
plt.xlabel('Characteristic Length (m)')
plt.ylabel('Pressure Error (Pa)')
plt.grid(which='both',axis='both')

# Add subplot for vertical displacement
fig.add_subplot(1,2,2)

# Plot displacement
# xData=[]
# yData=[]
for i in range(0,len(dt)):
	fileName=str(parentDirectory)+"/export/terzaghiErrorNorm_dt="+dt[i]+"_"+gridType+"-grid.txt"
	charLength=np.loadtxt(fname=fileName,usecols=(0))
	displacementError=np.loadtxt(fname=fileName,usecols=(2))
	plt.loglog(charLength,displacementError,'--o',color=colors[i],ms=5,mec='k',mew=0.5)
# 	data=np.loadtxt(fname=fileName)
# 	for i in range(0,pointsUsed):
# 		row=[]
# 		row.append(1)
# 		row.append(math.log10(data[i,0]))
# 		xData.append(row)
# 		yData.append(math.log10(data[i,2]))

# # Determine coefficients for linear regression
# xDataT=xData
# xData=np.transpose(xData)
# yData=np.transpose(yData)
# Coef=np.matmul(xData,xDataT)
# Coef=np.linalg.inv(Coef)
# yData=np.matmul(xData,yData)
# Coef=np.matmul(Coef,yData)

# # Plot linear regression
# hReg=[]
# linReg=[]
# for i in range(0,len(charLength)):
# 	hReg.append(charLength[i])
# 	linReg.append(lin(charLength[i]))
# plt.loglog(hReg,linReg,':k',linewidth=1.2,label='Lin.Reg.('+str(round(Coef[1],1))+')')

# Plot triangle
e1=1.05*displacementError[0]
h1=0.95*charLength[0]
h2=0.95*charLength[2]
e2=e1*((h2/h1)**2)
horCathetusX=[]
horCathetusY=[]
horCathetusX.append(h1)
horCathetusX.append(h2)
horCathetusY.append(e1)
horCathetusY.append(e1)
plt.loglog(horCathetusX,horCathetusY,'-k')
verCathetusX=[]
verCathetusY=[]
verCathetusX.append(h2)
verCathetusX.append(h2)
verCathetusY.append(e1)
verCathetusY.append(e2)
plt.loglog(verCathetusX,verCathetusY,'-k')
hypotenuseX=[]
hypotenuseY=[]
hypotenuseX.append(h1)
hypotenuseX.append(h2)
hypotenuseY.append(e1)
hypotenuseY.append(e2)
plt.loglog(hypotenuseX,hypotenuseY,'-k')

# Set axes' scale and limits
axes=plt.gca()
axes.set_xscale("log", nonposx='clip')
axes.set_yscale("log", nonposy='clip')
axes.text((h1+h2*2)/3,e1*1.05,"1")
axes.text(h2*0.9,(e1+e2*2)/3,"2")

# Set axes' labels
plt.xlabel('Characteristic Length (m)')
plt.ylabel('Vertical displacement error (m)')
plt.grid(which='both',axis='both')

# Add figure's legend
fig.legend(loc='upper center',ncol=4)

# Save figure
plt.savefig(plotName)

print("Plotted convergence results")