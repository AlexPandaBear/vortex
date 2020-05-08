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
X = [1., -1., 0., 0.]
Y = [0., 0., 1., -1.]
C = [1., 1., 1., 2.]
R = [0.01 for x in X]

#Time stepping parameters
t0 = 0.
tEnd = 100.
nb_steps = 100
nb_threads = 4

#Methods to compare
useEuler = True
useRK4 = True
useSV = True

trLen = 50 #Trace length in animation


#%% SCRIPT

nb_vtx = len(X)
H_tot = []
deltaAlpha = 1/trLen

if useEuler:
    sm_euler = vtx.SM()
    for i in range(nb_vtx):
        sm_euler.addVtx(X[i], Y[i], C[i], R[i])
    sm_euler.buildTimeSample(t0, tEnd, nb_steps)
    sm_euler.chooseNumericalMethod("euler")
    sm_euler.sim(nb_threads)
    T = sm_euler.getTimeList()
    H_euler = sm_euler.computeHamiltonianEvolution(nb_threads)
    H_tot += H_euler

if useRK4:
    sm_rk4 = vtx.SM()
    for i in range(nb_vtx):
        sm_rk4.addVtx(X[i], Y[i], C[i], R[i])
    sm_rk4.buildTimeSample(t0, tEnd, nb_steps)
    sm_rk4.chooseNumericalMethod("rk4")
    sm_rk4.sim(nb_threads)
    T = sm_rk4.getTimeList()
    H_rk4 = sm_rk4.computeHamiltonianEvolution(nb_threads)
    H_tot += H_rk4

if useSV:
    sm_sv = vtx.SM()
    for i in range(nb_vtx):
        sm_sv.addVtx(X[i], Y[i], C[i], R[i])
    sm_sv.buildTimeSample(t0, tEnd, nb_steps)
    sm_sv.chooseNumericalMethod("sv")
    sm_sv.sim(nb_threads)
    T = sm_sv.getTimeList()
    H_sv = sm_sv.computeHamiltonianEvolution(nb_threads)
    H_tot += H_sv

if useEuler or useRK4 or useSV:
    max = max(H_tot)
    min = min(H_tot)


plt.figure()

for t in range(nb_steps+1):
    plt.subplot(1,2,1)
    plt.cla()

    if useEuler:
        plt.plot(T, H_euler, label="euler", color='b')

    if useRK4:
        plt.plot(T, H_rk4, label="rk4", color='r')

    if useSV:
        plt.plot(T, H_sv, label="sv", color='g')

    if useEuler or useRK4 or useSV:
        plt.plot([T[t],T[t]], [min,max], color='k', linestyle='--')

    plt.xlabel("time (sec)")
    plt.ylabel("Hamiltonian")
    plt.title("Evolution of the Hamiltonian over the simulation")
    plt.legend()
    plt.show(block=False)

    plt.subplot(1,2,2)
    plt.show(block=False)

    if useEuler:
        X_euler = [sm_euler.getVtxXs(i) for i in range(nb_vtx)]
        Y_euler = [sm_euler.getVtxYs(i) for i in range(nb_vtx)]

    if useRK4:
        X_rk4 = [sm_rk4.getVtxXs(i) for i in range(nb_vtx)]
        Y_rk4 = [sm_rk4.getVtxYs(i) for i in range(nb_vtx)]

    if useSV:
        X_sv = [sm_sv.getVtxXs(i) for i in range(nb_vtx)]
        Y_sv = [sm_sv.getVtxYs(i) for i in range(nb_vtx)]

    b = t+1
    a = b - trLen

    if a < 0:
        a = 0

    plt.cla()
    plt.xlim([-2,2])
    plt.ylim([-2,2])
    plt.grid()

    for i in range(nb_vtx):
        if useEuler:
            plt.scatter(X_euler[i][t], Y_euler[i][t], color='b')
            X_plot = X_euler[i][a:b]
            Y_plot = Y_euler[i][a:b]
            for j in range(b-a):
                plt.plot(X_plot[j:j+2], Y_plot[j:j+2], color='b', alpha=(j+1)*deltaAlpha)

        if useRK4:
            plt.scatter(X_rk4[i][t], Y_rk4[i][t], color='r')
            X_plot = X_rk4[i][a:b]
            Y_plot = Y_rk4[i][a:b]
            for j in range(b-a):
                plt.plot(X_plot[j:j+2], Y_plot[j:j+2], color='r', alpha=(j+1)*deltaAlpha)

        if useSV:
            plt.scatter(X_sv[i][t], Y_sv[i][t], color='g')
            X_plot = X_sv[i][a:b]
            Y_plot = Y_sv[i][a:b]
            for j in range(b-a):
                plt.plot(X_plot[j:j+2], Y_plot[j:j+2], color='g', alpha=(j+1)*deltaAlpha)

    plt.draw()
    plt.pause(0.5)

plt.show()
