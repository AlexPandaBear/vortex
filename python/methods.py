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

x_max = max(X)
y_max = max(Y)

#Time stepping parameters
t0 = 0.
tEnd = 5000.
nb_steps = 4000
nb_threads = 4

#Methods to compare
useEuler = True
useRK4 = True
useEulerA = False
useEulerB = True
useSV = True
useSVI = False

methods = ['euler', 'rk4', 'eulerA', 'eulerB', 'sv', 'svi', 'test', 'test2', 'test3']
use_methods = [useEuler, useRK4, useEulerA, useEulerB, useSV, useSVI, False, False, False]
methods_colors = ['b', 'r', 'orange', 'g', 'hotpink', 'dodgerblue', 'r', 'k', 'b']
nb_methods = len(methods)

#Plotting parameters
trLen = 20 #Trace length in animation

normalize_H = True
normalize_T = True

T_period_th = 8 * np.pi**2 # 2 pi R / [Gamma / (2 pi D)] = 8 pi^2 R^2 / Gamma


#%% SCRIPT

nb_vtx = len(X)
H_tot = []
deltaAlpha = 1/trLen

print("Starting simulation(s)...")

SM = []
H = []

for i in range(nb_methods):
    if use_methods[i]:
        sm = vtx.SM()
        for j in range(nb_vtx):
            sm.addVtx(X[j], Y[j], C[j], R[j], 0)
        sm.buildTimeSample(t0, tEnd, nb_steps)
        sm.chooseNumericalMethod(methods[i])
        sm.sim(nb_threads)

        T = sm.getTimeList()
        if normalize_T:
            T = [t/T_period_th for t in T]

        H_tmp = sm.computeHamiltonianEvolution(nb_threads)
        if normalize_H:
            H_tmp = [h/H_tmp[0] for h in H_tmp]

        SM.append(sm)
        H.append(H_tmp)
        H_tot += H_tmp

    else:
        SM.append(False)
        H.append(False)

print("Simulation(s) teminated")

h_max = max(H_tot)
h_min = min(H_tot)

print("Plotting...")

fig, (ax1, ax2) = plt.subplots(1,2)
X = []
Y = []
X_th = [np.cos(theta) for theta in np.linspace(0, 2*np.pi, 1000)]
Y_th = [np.sin(theta) for theta in np.linspace(0, 2*np.pi, 1000)]

for t in range(nb_steps+1):
    ax1.cla()

    for i in range(nb_methods):
        if use_methods[i]:
            ax1.plot(T,H[i], label=methods[i], color=methods_colors[i])

    ax1.plot([T[t],T[t]], [h_min,h_max], color='k', linestyle='--')

    if normalize_T:
        ax1.set_xlabel(r'Normalized time $t^* = t \; / \; T_{per}$')
    else:
        ax1.set_xlabel("time (sec)")

    if normalize_H:
        ax1.set_ylabel(r'Normalized Hamiltonian $H^* = H \; / \; H_0$')
    else:
        ax1.set_ylabel("Hamiltonian")

    ax1.set_title("Evolution of the Hamiltonian over the simulation")
    ax1.legend()
    plt.show(block=False)

    for i in range(nb_methods):
        if use_methods[i]:
            X.append([SM[i].getVtxXs(j) for j in range(nb_vtx)])
            Y.append([SM[i].getVtxYs(j) for j in range(nb_vtx)])
        else:
            X.append(False)
            Y.append(False)

    b = t+1
    a = b - trLen

    if a < 0:
        a = 0

    ax2.cla()

    local_x_max = 0
    local_y_max = 0

    for i in range(nb_vtx):
        for k in range(len(methods)):
            if use_methods[k]:
                local_x_max = max([local_x_max, X[k][i][t]])
                local_y_max = max([local_y_max, Y[k][i][t]])

    x_max = max([x_max, local_x_max])
    y_max = max([y_max, local_y_max])

    lim = 1.1 * max([x_max, y_max])
    ax2.set_xlim([-lim,lim])
    ax2.set_ylim([-lim,lim])

    ax2.grid()
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")

    ax2.plot(X_th, Y_th, color='k', linestyle='--')

    for i in range(nb_vtx):
        for k in range(nb_methods):
            if use_methods[k]:
                ax2.scatter(X[k][i][t], Y[k][i][t], color=methods_colors[k])
                X_plot = X[k][i][a:b]
                Y_plot = Y[k][i][a:b]
                for j in range(b-a):
                    ax2.plot(X_plot[j:j+2], Y_plot[j:j+2], color=methods_colors[k], alpha=(j+1)*deltaAlpha)

    plt.draw()
    plt.pause(0.5)

plt.show()

print("--------------------------------------------------------END-")
