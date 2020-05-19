print("------------------------------------------------------------")
print("------------------   VORTEX METHOD CODE   ------------------")
print("--------------------   Method Testing   --------------------")
print("-------   Alexandre DUTKA - ISAE-SUPAERO - 05/2020   -------")
print("------------------------------------------------------------")


#%% IMPORTS

import numpy as np
import _vortex as vtx
import matplotlib.pyplot as plt


#%% PARAMETERS

#Vortex definition
X = [1., -1.]
Y = [0., 0.]
C = [1., 1.]
R = [0.01 for x in X]

#Time stepping parameters
t0 = 0.
tEnd = 10000.
nb_steps = [10000, 20000, 40000, 80000, 160000]
nb_threads = 4
method = "euler"

trLen = 20 #Trace length in animation


#%% SCRIPT

nb_vtx = len(X)

print("Starting simulation(s)...")

SM = []
for nbs in nb_steps:
    sm = vtx.SM()
    for i in range(nb_vtx):
        sm.addVtx(X[i], Y[i], C[i], R[i], 0)
    sm.buildTimeSample(t0, tEnd, nbs)
    sm.chooseNumericalMethod(method)
    sm.sim(nb_threads)
    SM.append(sm)

print("Simulation(s) teminated")
print("Plotting...")

plt.figure()

for sm in SM:
    T = sm.getTimeList()
    H = sm.computeHamiltonianEvolution(nb_threads)
    plt.plot(T, H, label="{} steps".format(len(T)-1))

plt.plot([T[0], T[-1]], [H[0], H[0]], color='k', linestyle='--')
plt.xlabel("time (sec)")
plt.ylabel("Hamiltonian")
plt.title("Evolution of the Hamiltonian over the simulation")
plt.legend()
plt.show()

print("--------------------------------------------------------END-")
