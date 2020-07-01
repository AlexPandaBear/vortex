print("------------------------------------------------------------")
print("------------------   VORTEX METHOD CODE   ------------------")
print("-------------   KELVIN-HELMHOLTZ INSTABILITY   -------------")
print("----------------------   Simulation   ----------------------")
print("-------   Alexandre DUTKA - ISAE-SUPAERO - 06/2020   -------")
print("------------------------------------------------------------")


#%% IMPORTS

print("Importing libraries")

import numpy as np
import _vortex as vtx
import matplotlib.pyplot as plt
import time


#%% PARAMETERS

print("Loading parameters")

width = 15. #size (x-axis) of the "main box"
height = 10. #size (y-axis) of the "main box"

U0 = 1. #initial velocity of the top layer (bottom one's is -U0)

RSLT = 0.2 #Relative Shear Layer Thickness (normalized by the height parameter)
RIHPA = 0.02 #Relative Initial Harmonic Perturbation Amplitude (also normalized by height)

nx = 50 #number of columns for the vtx matrix (initialization)
ny = 50 #number of lines for the vtx matrix (initialization)

#integration parameters
t0 = 0.
tEnd = 20.
nb_steps = 1000

#Numerical methods implemented : Explicit Euler (euler), Runge-Kutta 4 (rk4), Asymetrical Euler-A (eulerA), Asymetrical Euler-B (eulerB), Stormer-Verlet (sv) and Stormer-Verlet Inverse (svi)
temporalIntegrationMethod = "rk4"

periodicity = True #has to be True to properly compute the Kelvin-Helmholtz instability

dx = width/(nx-1) #horizontal space between two vtx at t0
dy = height/(ny-1) #vertical space between two vtx at t0

SLT = 0.5*height*RSLT #Shear Layer Thickness
IHPA = 0.5*height*RIHPA #Initial Harmonic Perturbation Amplitude

regRadius = max(dx, dy) #Regularization radius of the vortices


#%% INSTRUCTIONS

print("Reading instructions")

ignoreNullCircVtx = True #if True, vortices with null circulation in the initialization matix will be deleted
showRealVtxNumber = True #if True, prints the number of vortices really used for computation

plotInitialProfiles = True #if True, plots the velocity and vorticity profiles before stating computation

showComputationTime = True #if True, prints the computation time

saveSim = True #if False, once this script terminates all the simulation will be lost
dataFolder = "../data" #folder where simulation files will be stored (default is "../data")
saveFile = "rk4" #name of the simulation files (if it matches an existing file, it will be overwritten)

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
		yWithPert = y + IHPA * np.sin(2*np.pi * x/width) * np.cos(np.pi * y/height)
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
