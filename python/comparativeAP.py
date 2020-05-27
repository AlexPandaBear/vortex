print("------------------------------------------------------------")
print("------------------   VORTEX METHOD CODE   ------------------")
print("-------------   Comparative After-Processing   -------------")
print("-------   Alexandre DUTKA - ISAE-SUPAERO - 05/2020   -------")
print("------------------------------------------------------------")


#%% IMPORTS

print("Importing libraries")

import matplotlib.pyplot as plt
import time
import math
import numpy as np
import _vortex as vtx


#%% INSTRUCTIONS

print("Reading instructions")

loadNewData = True #if False, previously loaded data will be used
dataFolder = "../data"
dataFiles = ["test", "test_rk4", "test_sv"]
simNames = ["Euler", "Runge-Kutta 4", "Stormer-Verlet"]

compareParameters = True

compareVtxConfigs = False
colors_vtxConfigs = ["dodgerblue", "hotpink"]
vtxConfigs_reframe = True
steps_vtxConfigs = [i for i in range(1001) if i%50==0]

compareHamiltonians = True

numberOfThreads = 4


#%% FUNCTIONS

print("Loading functions")

#Loads the necessary workspace, saved in the file file.wsdata (associated to the simulation file.simdata)
def loadWorkspace(file, WS):
    f = open(file, "r")
    ws = []
    width = float(f.readline())
    ws.append(width)
    height = float(f.readline())
    ws.append(height)
    ws.append(float(f.readline()))
    RSLT = float(f.readline())
    ws.append(RSLT)
    ws.append(float(f.readline()))
    nx = int(f.readline())
    ws.append(nx)
    ny = int(f.readline())
    ws.append(ny)
    ws.append(float(f.readline()))
    ws.append(float(f.readline()))
    ws.append(int(f.readline()))
    ws.append(f.readline().rstrip())
    periodicity = f.readline().rstrip()
    if periodicity == "True":
        ws.append(True)
    else:
        ws.append(False)
    ws.append(width/(nx-1))
    ws.append(height/(ny-1))
    ws.append(0.5*height*RSLT)
    ws.append(float(f.readline()))
    f.close()
    WS.append(ws)

def showSimParameters(W):
    print("------------------------------------------------------------")
    firstLine = "Simulation  "
    for name in simNames:
        firstLine += " | " + name + " "*(20-len(name))
    print(firstLine)
    print("------------------------------------------------------------")
    lines = ["width", "height", "U0", "RSLT", "RIHPA", "nx", "ny", "t0", "tEnd", "nb_steps", "method", "periodicity", "dx", "dy", "SLT", "regRadius"]
    for i in range(len(lines)):
        l = lines[i] + " "*(12-len(lines[i]))
        for w in W:
            l += " | " + str(w[i]) + " "*(20-len(str(w[i])))
        print(l)
    print("------------------------------------------------------------")


#Plots the vortices' positions at step step for each simulation
def compareVtxPositions(S, W, step):
    plt.figure()

    for k in range(len(S)):
        s = S[k]
        X = s.getXsAt(step)
        Y = s.getYsAt(step)
        C = s.getCirculationsAt(step)
        nb_vtx = len(X)

        if vtxConfigs_reframe:
            width = WS[k][0]

            for i in range(nb_vtx):
            	if X[i] > width:
            		X[i] -= width
            	elif X[i] < 0.:
            		X[i] += width

        circMax = max([abs(circ) for circ in C])

        nameNotUsed = True
        for i in range(nb_vtx):
            if nameNotUsed and abs(C[i]) == circMax:
                plt.scatter(X[i], Y[i], 100*abs(C[i])/circMax, c=colors_vtxConfigs[k], label=simNames[k])
                nameNotUsed = False
            else:
                plt.scatter(X[i], Y[i], 100*abs(C[i])/circMax, c=colors_vtxConfigs[k])


    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.legend()
    plt.title("Vortex configurations at step {}".format(step))
    plt.show()

#Plots the evolution of the Hamiltonian of the system during each simulation stored in the sim. manager list S
def compareHamiltonianEvolutions(S, nbThreads):
    plt.figure()

    for k in range(len(S)):
        s = S[k]
        T = s.getTimeList()
        H = s.computeHamiltonianEvolution(nbThreads)
        plt.plot([T[0], T[-1]], [H[0], H[0]], c='k', label="H0 "+simNames[k], linestyle='--')
        plt.plot(T, H, label=simNames[k])

    plt.xlabel("time (sec)")
    plt.ylabel("Hamiltonian")
    plt.legend()
    plt.title("Evolution of the Hamiltonian over the simulation")
    plt.show()


#%% SCRIPT

if loadNewData:
    SM = []
    WS = []

    for file in dataFiles:
        print("Loading simulation data from file " + dataFolder + "/" + file + ".simdata")
        sm = vtx.SM()
        sm.loadSim(dataFolder + "/" + file + ".simdata")
        sm.setName(file)
        SM.append(sm)

        print("Loading workspace data from file " + dataFolder + "/" + file + ".wsdata")
        loadWorkspace(dataFolder + "/" + file + ".wsdata", WS)

    print("Data loaded successfully")

if compareParameters:
    print("Simulation parameters")
    showSimParameters(WS)

if compareVtxConfigs:
    for step in steps_vtxConfigs:
        print("Plotting vortex positions comparaison at step {}".format(step))
        compareVtxPositions(SM, WS, step)

if compareHamiltonians:
    print("Plotting Hamiltonian evolutions comparaison")
    compareHamiltonianEvolutions(SM, numberOfThreads)

print("--------------------------------------------------------END-")
