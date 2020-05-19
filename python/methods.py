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
nb_steps = 10000
nb_threads = 4

#Methods to compare
useEuler = True
useRK4 = True
useEulerA = True
useEulerB = False
useSV = True
useSVI = True

trLen = 20 #Trace length in animation


#%% SCRIPT

nb_vtx = len(X)
H_tot = []
deltaAlpha = 1/trLen

print("Starting simulation(s)...")

if useEuler:
    sm_euler = vtx.SM()
    for i in range(nb_vtx):
        sm_euler.addVtx(X[i], Y[i], C[i], R[i], 0)
    sm_euler.buildTimeSample(t0, tEnd, nb_steps)
    sm_euler.chooseNumericalMethod("euler")
    sm_euler.sim(nb_threads)
    T = sm_euler.getTimeList()
    H_euler = sm_euler.computeHamiltonianEvolution(nb_threads)
    H_tot += H_euler

if useRK4:
    sm_rk4 = vtx.SM()
    for i in range(nb_vtx):
        sm_rk4.addVtx(X[i], Y[i], C[i], R[i], 0)
    sm_rk4.buildTimeSample(t0, tEnd, nb_steps)
    sm_rk4.chooseNumericalMethod("rk4")
    sm_rk4.sim(nb_threads)
    T = sm_rk4.getTimeList()
    H_rk4 = sm_rk4.computeHamiltonianEvolution(nb_threads)
    H_tot += H_rk4

if useEulerA:
    sm_eulerA = vtx.SM()
    for i in range(nb_vtx):
        sm_eulerA.addVtx(X[i], Y[i], C[i], R[i], 0)
    sm_eulerA.buildTimeSample(t0, tEnd, nb_steps)
    sm_eulerA.chooseNumericalMethod("eulerA")
    sm_eulerA.sim(nb_threads)
    T = sm_eulerA.getTimeList()
    H_eulerA = sm_eulerA.computeHamiltonianEvolution(nb_threads)
    H_tot += H_eulerA

if useEulerB:
    sm_eulerB = vtx.SM()
    for i in range(nb_vtx):
        sm_eulerB.addVtx(X[i], Y[i], C[i], R[i], 0)
    sm_eulerB.buildTimeSample(t0, tEnd, nb_steps)
    sm_eulerB.chooseNumericalMethod("eulerB")
    sm_eulerB.sim(nb_threads)
    T = sm_eulerB.getTimeList()
    H_eulerB = sm_eulerB.computeHamiltonianEvolution(nb_threads)
    H_tot += H_eulerB

if useSV:
    sm_sv = vtx.SM()
    for i in range(nb_vtx):
        sm_sv.addVtx(X[i], Y[i], C[i], R[i], 0)
    sm_sv.buildTimeSample(t0, tEnd, nb_steps)
    sm_sv.chooseNumericalMethod("sv")
    sm_sv.sim(nb_threads)
    T = sm_sv.getTimeList()
    H_sv = sm_sv.computeHamiltonianEvolution(nb_threads)
    H_tot += H_sv

if useSVI:
    sm_svi = vtx.SM()
    for i in range(nb_vtx):
        sm_svi.addVtx(X[i], Y[i], C[i], R[i], 0)
    sm_svi.buildTimeSample(t0, tEnd, nb_steps)
    sm_svi.chooseNumericalMethod("sv")
    sm_svi.sim(nb_threads)
    T = sm_svi.getTimeList()
    H_svi = sm_svi.computeHamiltonianEvolution(nb_threads)
    H_tot += H_svi

print("Simulation(s) teminated")

if useEuler or useRK4 or useEulerA or useEulerB or useSV or useVI:
    max = max(H_tot)
    min = min(H_tot)

print("Plotting...")

fig, (ax1, ax2) = plt.subplots(1,2)

for t in range(nb_steps+1):
    ax1.cla()

    if useEuler:
        ax1.plot(T, H_euler, label="euler", color='b')

    if useRK4:
        ax1.plot(T, H_rk4, label="rk4", color='r')

    if useEulerA:
        ax1.plot(T, H_eulerA, label="eulerA", color='orange')

    if useEulerB:
        ax1.plot(T, H_eulerB, label="eulerB", color='violet')

    if useSV:
        ax1.plot(T, H_sv, label="sv", color='g')

    if useSVI:
        ax1.plot(T, H_svi, label="svi", color='c')

    if useEuler or useRK4 or useSV or useSVI:
        ax1.plot([T[t],T[t]], [min,max], color='k', linestyle='--')

    ax1.set_xlabel("time (sec)")
    ax1.set_ylabel("Hamiltonian")
    ax1.set_title("Evolution of the Hamiltonian over the simulation")
    ax1.legend()
    plt.show(block=False)

    if useEuler:
        X_euler = [sm_euler.getVtxXs(i) for i in range(nb_vtx)]
        Y_euler = [sm_euler.getVtxYs(i) for i in range(nb_vtx)]

    if useRK4:
        X_rk4 = [sm_rk4.getVtxXs(i) for i in range(nb_vtx)]
        Y_rk4 = [sm_rk4.getVtxYs(i) for i in range(nb_vtx)]

    if useEulerA:
        X_eulerA = [sm_eulerA.getVtxXs(i) for i in range(nb_vtx)]
        Y_eulerA = [sm_eulerA.getVtxYs(i) for i in range(nb_vtx)]

    if useEulerB:
        X_eulerB = [sm_eulerB.getVtxXs(i) for i in range(nb_vtx)]
        Y_eulerB = [sm_eulerB.getVtxYs(i) for i in range(nb_vtx)]

    if useSV:
        X_sv = [sm_sv.getVtxXs(i) for i in range(nb_vtx)]
        Y_sv = [sm_sv.getVtxYs(i) for i in range(nb_vtx)]

    if useSVI:
        X_svi = [sm_svi.getVtxXs(i) for i in range(nb_vtx)]
        Y_svi = [sm_svi.getVtxYs(i) for i in range(nb_vtx)]

    b = t+1
    a = b - trLen

    if a < 0:
        a = 0

    ax2.cla()
    ax2.set_xlim([-2,2])
    ax2.set_ylim([-2,2])
    ax2.grid()
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")

    for i in range(nb_vtx):
        if useEuler:
            ax2.scatter(X_euler[i][t], Y_euler[i][t], color='b')
            X_plot = X_euler[i][a:b]
            Y_plot = Y_euler[i][a:b]
            for j in range(b-a):
                ax2.plot(X_plot[j:j+2], Y_plot[j:j+2], color='b', alpha=(j+1)*deltaAlpha)

        if useRK4:
            ax2.scatter(X_rk4[i][t], Y_rk4[i][t], color='r')
            X_plot = X_rk4[i][a:b]
            Y_plot = Y_rk4[i][a:b]
            for j in range(b-a):
                ax2.plot(X_plot[j:j+2], Y_plot[j:j+2], color='r', alpha=(j+1)*deltaAlpha)

        if useEulerA:
            ax2.scatter(X_eulerA[i][t], Y_eulerA[i][t], color='orange')
            X_plot = X_eulerA[i][a:b]
            Y_plot = Y_eulerA[i][a:b]
            for j in range(b-a):
                ax2.plot(X_plot[j:j+2], Y_plot[j:j+2], color='orange', alpha=(j+1)*deltaAlpha)

        if useEulerB:
            ax2.scatter(X_eulerB[i][t], Y_eulerB[i][t], color='violet')
            X_plot = X_eulerB[i][a:b]
            Y_plot = Y_eulerB[i][a:b]
            for j in range(b-a):
                ax2.plot(X_plot[j:j+2], Y_plot[j:j+2], color='violet', alpha=(j+1)*deltaAlpha)

        if useSV:
            ax2.scatter(X_sv[i][t], Y_sv[i][t], color='g')
            X_plot = X_sv[i][a:b]
            Y_plot = Y_sv[i][a:b]
            for j in range(b-a):
                ax2.plot(X_plot[j:j+2], Y_plot[j:j+2], color='g', alpha=(j+1)*deltaAlpha)

        if useSVI:
            ax2.scatter(X_svi[i][t], Y_svi[i][t], color='c')
            X_plot = X_svi[i][a:b]
            Y_plot = Y_svi[i][a:b]
            for j in range(b-a):
                ax2.plot(X_plot[j:j+2], Y_plot[j:j+2], color='c', alpha=(j+1)*deltaAlpha)

    plt.draw()
    plt.pause(0.5)

plt.show()

print("--------------------------------------------------------END-")
