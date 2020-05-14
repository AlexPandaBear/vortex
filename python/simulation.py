print("------------------------------------------------------------")
print("------------------   VORTEX METHOD CODE   ------------------")
print("----------------------   Simulation   ----------------------")
print("-------   Alexandre DUTKA - ISAE-SUPAERO - 05/2020   -------")
print("------------------------------------------------------------")


#%% IMPORTS

print("Importing libraries")

import numpy as np
import _vortex as vtx
import matplotlib.pyplot as plt
import time


#%% PARAMETERS

print("Loading parameters")

width = 10.
height = 10.

U0 = 1. #initial velocity of the top layer (bottom one's is -U0)

RSLT = 0.3 #relative shear layer thickness
RIHPA = 0.3 #relative initial harmonic perturbation amplitude

nx = 50 #number of columns for the vtx matrix
ny = 50 #number of lines for the vtx matrix

t0 = 0.
tEnd = 40.
nb_steps = 200

#Numerical methods implemented : Explicit Euler (euler), Runge-Kutta 4 (rk4) and Stormer-Verlet (sv)
temporalIntegrationMethod = "rk4"

periodicity = True #computes mvt considering the periodicity of the flow (4 neighbour cells used)

dx = width/(nx-1) #horizontal space between two vtx at t0
dy = height/(ny-1) #vertical space between two vtx at t0

SLT = 0.5*height*RSLT #shear layer thickness

regRadius = max(dx, dy)


#%% INSTRUCTIONS

print("Reading instructions")

ignoreNullCircVtx = True
showRealVtxNumber = True

plotInitialProfiles = True

showComputationTime = True

saveSim = True
dataFolder = "../data"
saveFile = "test3"

numberOfThreads = 4


#%% FUNCTIONS

print("Loading functions")

#Computes velocity along x axis given the y coordinate and the shear layer thickness e
def u(y, e):
	if y > e:
		return U0
	elif y < -e:
		return -U0
	else:
		return U0*np.sin(0.5*np.pi*y/e)

#Computes circulation of a vortex given the x and y coordinates and the shear layer thickness e
def computeVtxCircAt(x, y, e):
	"""
	omega = - (u(y+dy, e) - u(y-dy, e))/(2*dy)
	gamma = omega*dx*dy
	return gamma
	"""
	return -0.5*(u(y+dy, e) - u(y-dy, e))*dx

#Plots the velocity and vorticity profiles at the beginning of the simulation s (ie t = t0)
def initialProfiles():
	plt.figure()

	plt.subplot(1,2,1)
	Y = [-height/2. + i*dy for i in range(ny)]
	U = [u(y, SLT) for y in Y]
	plt.plot([0,0],[-0.5*height, 0.5*height], 'k', linewidth=1, linestyle='--')
	plt.plot([-U0,U0],[0,0], 'k', linewidth=1, linestyle='--')
	plt.plot(U, Y)
	plt.xlabel("U (m/s)")
	plt.ylabel("y (m)")

	plt.subplot(1,2,2)
	O = [U[i+1] - U[i-1] for i in range(1,ny-1)]
	O.insert(0,0.)
	O.append(0.)
	plt.plot([0,0],[-0.5*height, 0.5*height], 'k', linewidth=1, linestyle='--')
	plt.plot(O, Y)
	plt.xlabel(r"$\omega$ (rad/s)")

	plt.show()


#Saves the workspace associated to the simulation in the file <file>.wsdata
def saveWorkspace(file):
    f = open(file, "w")
    f.writelines(str(width) + "\n")
    f.writelines(str(height) + "\n")
    f.writelines(str(U0) + "\n")
    f.writelines(str(RSLT) + "\n")
    f.writelines(str(RIHPA) + "\n")
    f.writelines(str(nx) + "\n")
    f.writelines(str(ny) + "\n")
    f.writelines(str(t0) + "\n")
    f.writelines(str(tEnd) + "\n")
    f.writelines(str(nb_steps) + "\n")
    f.writelines(temporalIntegrationMethod + "\n")
    f.writelines(str(periodicity) + "\n")
    f.writelines(str(regRadius) + "\n")
    f.close()


#%% SCRIPT

print("Creating a new simulation environment")
sm = vtx.SM()

print("Intializing simulation")
for i in range(nx):
	for j in range(ny):
		x = width * i/(nx)
		y = height * (0.5 - j/(ny-1))
		circ = computeVtxCircAt(x, y, SLT)
		yWithPert = y + RIHPA * np.sin(2*np.pi * x/width) * np.cos(np.pi * y/height)
		fluid = 0
		if y < 0.:
			fluid = 1
		if ignoreNullCircVtx:
			if circ!=0.:
				sm.addVtx(x, yWithPert, circ, regRadius, fluid)
		else:
			sm.addVtx(x, yWithPert, circ, regRadius, fluid)

sm.buildTimeSample(t0, tEnd, nb_steps)
sm.setXPeriodicityTo(periodicity, width)
sm.chooseNumericalMethod(temporalIntegrationMethod)

if plotInitialProfiles:
    print("Plotting initial velocity and vorticity profiles")
    initialProfiles()

print("Starting simulation...")
tStart = time.time()
sm.sim(numberOfThreads)
print("Simulation finished successfully")

if showRealVtxNumber:
	print("Number of vortices used for computation : {}".format(sm.getNbVtx()))

if showComputationTime:
	print("Computation time : {} sec".format(time.time() - tStart))

if saveSim:
    if saveFile == "":
        print("Unable to save data (bad saveFile name)")
    else:
        print("Saving simulation data in file " + dataFolder + "/" + saveFile + ".simdata")
        sm.saveSim(dataFolder + "/" + saveFile + ".simdata")
        print("Saving workspace data in file " + dataFolder + "/" + saveFile + ".wsdata")
        saveWorkspace(dataFolder + "/" + saveFile + ".wsdata")
        print("Data saved successfully")

print("--------------------------------------------------------END-")
