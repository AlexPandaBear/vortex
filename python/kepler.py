print("------------------------------------------------------------")
print("----------------   2D KEPLER PROBLEM CODE   ----------------")
print("-------   Alexandre DUTKA - ISAE-SUPAERO - 06/2020   -------")
print("------------------------------------------------------------")

import numpy as np
import matplotlib.pyplot as plt


#%% CONSTANTS

G = 6.674e-11
T_year = 365.25*86400.
AU = 150000000000. #Sun-Earth distance


#%% PARAMETERS

#initial conditions - positions
X0 = [0., AU]
Y0 = [0., 0.]

#initial conditions - velocities
U0 = [0., 0.]
V0 = [0., 30000.]

MASS = [1.989e30, 5.972e24]
NAME = ['Sun','Earth']

#integration parameters
t0 = 0.
tEnd = 3.*T_year
nb_steps = 300

#methods to compare
useEuler = True
useEulerA = True
useEulerB = True
useSV = True
useSVI = False
useRK4 = True


#%% SCRIPT

methods = ["Explicit Euler", "EulerA", "EulerB", "Stormer-Verlet", "Inverse Stormer-Verlet", "Runge-Kutta 4"]
useMethods = [useEuler, useEulerA, useEulerB, useSV, useSVI, useRK4]
nb_methods = len(useMethods)
methodColors = ['b', 'orange', 'hotpink', 'g', 'dodgerblue', 'r']

trLen = 100

nb_planets = len(X0)
dt = (tEnd - t0)/nb_steps

Xs = [np.zeros((nb_planets, nb_steps+1)) for i in range(nb_methods)]
Ys = [np.zeros((nb_planets, nb_steps+1)) for i in range(nb_methods)]
Us = [np.zeros((nb_planets, nb_steps+1)) for i in range(nb_methods)]
Vs = [np.zeros((nb_planets, nb_steps+1)) for i in range(nb_methods)]

for m in range(nb_methods):
    Xs[m][:,0] = X0
    Ys[m][:,0] = Y0
    Us[m][:,0] = U0
    Vs[m][:,0] = V0

def Fx(method, planet, step):
    f = 0.
    for i in range(nb_planets):
        if i != planet:
            dist3 = ((Xs[method][i,step] - Xs[method][planet,step])**2 + (Ys[method][i,step] - Ys[method][planet,step])**2)**1.5
            f += MASS[i]*(Xs[method][i,step] - Xs[method][planet,step])/dist3
    f *= G*MASS[planet]
    return f

def Fy(method, planet, step):
    f = 0.
    for i in range(nb_planets):
        if i != planet:
            dist3 = ((Xs[method][i,step] - Xs[method][planet,step])**2 + (Ys[method][i,step] - Ys[method][planet,step])**2)**1.5
            f += MASS[i]*(Ys[method][i,step] - Ys[method][planet,step])/dist3
    f *= G*MASS[planet]
    return f

def Fx_SVI(planet, X_tmp, Y_tmp):
    f = 0.
    for i in range(nb_planets):
        if i != planet:
            dist3 = ((X_tmp[i] - X_tmp[planet])**2 + (Y_tmp[i] - Y_tmp[planet])**2)**1.5
            f += MASS[i]*(X_tmp[i] - X_tmp[planet])/dist3
    f *= G*MASS[planet]
    return f

def Fy_SVI(planet, X_tmp, Y_tmp):
    f = 0.
    for i in range(nb_planets):
        if i != planet:
            dist3 = ((X_tmp[i] - X_tmp[planet])**2 + (Y_tmp[i] - Y_tmp[planet])**2)**1.5
            f += MASS[i]*(Y_tmp[i] - Y_tmp[planet])/dist3
    f *= G*MASS[planet]
    return f

def computeEulerStep(step):
    for i in range(nb_planets):
        Xs[0][i,step+1] = Xs[0][i,step] + dt*Us[0][i,step]
        Ys[0][i,step+1] = Ys[0][i,step] + dt*Vs[0][i,step]
        Us[0][i,step+1] = Us[0][i,step] + dt*Fx(0, i, step)/MASS[i]
        Vs[0][i,step+1] = Vs[0][i,step] + dt*Fy(0, i, step)/MASS[i]

def computeEulerAStep(step):
    for i in range(nb_planets):
        Xs[1][i,step+1] = Xs[1][i,step] + dt*Us[1][i,step]
        Ys[1][i,step+1] = Ys[1][i,step] + dt*Vs[1][i,step]
    for i in range(nb_planets):
        Us[1][i,step+1] = Us[1][i,step] + dt*Fx(1, i, step+1)/MASS[i]
        Vs[1][i,step+1] = Vs[1][i,step] + dt*Fy(1, i, step+1)/MASS[i]

def computeEulerBStep(step):
    for i in range(nb_planets):
        Us[2][i,step+1] = Us[2][i,step] + dt*Fx(2, i, step)/MASS[i]
        Vs[2][i,step+1] = Vs[2][i,step] + dt*Fy(2, i, step)/MASS[i]
    for i in range(nb_planets):
        Xs[2][i,step+1] = Xs[2][i,step] + dt*Us[2][i,step+1]
        Ys[2][i,step+1] = Ys[2][i,step] + dt*Vs[2][i,step+1]

def computeSVStep(step):
    U_tmp = np.zeros(nb_planets)
    V_tmp = np.zeros(nb_planets)
    for i in range(nb_planets):
        U_tmp[i] = Us[3][i,step] + 0.5*dt*Fx(3, i, step)/MASS[i]
        V_tmp[i] = Vs[3][i,step] + 0.5*dt*Fy(3, i, step)/MASS[i]
    for i in range(nb_planets):
        Xs[3][i,step+1] = Xs[3][i,step] + dt*U_tmp[i]
        Ys[3][i,step+1] = Ys[3][i,step] + dt*V_tmp[i]
    for i in range(nb_planets):
        Us[3][i,step+1] = U_tmp[i] + 0.5*dt*Fx(3, i, step+1)/MASS[i]
        Vs[3][i,step+1] = V_tmp[i] + 0.5*dt*Fy(3, i, step+1)/MASS[i]

def computeSVIStep(step):
    X_tmp = np.zeros(nb_planets)
    Y_tmp = np.zeros(nb_planets)
    for i in range(nb_planets):
        X_tmp[i] = Xs[4][i,step] + 0.5*dt*Us[4][i,step]
        Y_tmp[i] = Ys[4][i,step] + 0.5*dt*Vs[4][i,step]
    for i in range(nb_planets):
        Us[4][i,step+1] = Us[4][i,step] + dt*Fx_SVI(i,X_tmp,Y_tmp)
        Vs[4][i,step+1] = Vs[4][i,step] + dt*Fy_SVI(i,X_tmp,Y_tmp)
    for i in range(nb_planets):
        Xs[4][i,step+1] = X_tmp[i] + 0.5*dt*Us[4][i,step+1]
        Ys[4][i,step+1] = Y_tmp[i] + 0.5*dt*Vs[4][i,step+1]

def computeRK4Step(step):
    X_tmp = np.zeros(nb_planets)
    Y_tmp = np.zeros(nb_planets)
    k1 = np.zeros((nb_planets,2))
    k2 = np.zeros((nb_planets,2))
    k3 = np.zeros((nb_planets,2))
    k4 = np.zeros((nb_planets,2))

    for i in range(nb_planets):
        k1[i,0] = Fx(5, i, step)/MASS[i]
        k1[i,1] = Fy(5, i, step)/MASS[i]

    for i in range(nb_planets):
        X_tmp[i] = Xs[5][i,step] + 0.5*dt*Us[5][i,step]
        Y_tmp[i] = Ys[5][i,step] + 0.5*dt*Vs[5][i,step]

    for i in range(nb_planets):
        k2[i,0] = Fx_SVI(i, X_tmp, Y_tmp)/MASS[i]
        k2[i,1] = Fy_SVI(i, X_tmp, Y_tmp)/MASS[i]

    for i in range(nb_planets):
        X_tmp[i] += 0.25*(dt**2)*k1[i,0]
        Y_tmp[i] += 0.25*(dt**2)*k1[i,1]

    for i in range(nb_planets):
        k3[i,0] = Fx_SVI(i, X_tmp, Y_tmp)/MASS[i]
        k3[i,1] = Fy_SVI(i, X_tmp, Y_tmp)/MASS[i]

    for i in range(nb_planets):
        X_tmp[i] = Xs[5][i,step] + 0.5*dt*Us[5][i,step] + 0.5*(dt**2)*k2[i,0]
        Y_tmp[i] = Ys[5][i,step] + 0.5*dt*Vs[5][i,step] + 0.5*(dt**2)*k2[i,1]

    for i in range(nb_planets):
        k4[i,0] = Fx_SVI(i, X_tmp, Y_tmp)/MASS[i]
        k4[i,1] = Fy_SVI(i, X_tmp, Y_tmp)/MASS[i]

    for i in range(nb_planets):
        Xs[5][i,step+1] = Xs[5][i,step] + dt*Us[5][i,step] + (dt**2)*(k1[i,0] + k2[i,0] + k3[i,0])/6.
        Ys[5][i,step+1] = Ys[5][i,step] + dt*Vs[5][i,step] + (dt**2)*(k1[i,1] + k2[i,1] + k3[i,1])/6.
        Us[5][i,step+1] = Us[5][i,step] + dt*(k1[i,0] + 2.*k2[i,0] + 2.*k3[i,0] + k4[i,0])/6.
        Vs[5][i,step+1] = Vs[5][i,step] + dt*(k1[i,1] + 2.*k2[i,1] + 2.*k3[i,1] + k4[i,1])/6.

for step in range(nb_steps):
    if useEuler:
        computeEulerStep(step)
    if useEulerA:
        computeEulerAStep(step)
    if useEulerB:
        computeEulerBStep(step)
    if useSV:
        computeSVStep(step)
    if useSVI:
        computeSVIStep(step)
    if useRK4:
        computeRK4Step(step)

def computeHamiltonianAt(method, step):
    T = 0.
    V = 0.
    for i in range(nb_planets):
        T += MASS[i]*(Us[method][i,step]**2 + Vs[method][i,step]**2)
        for j in range(i+1,nb_planets):
            dist = np.sqrt((Xs[method][i,step] - Xs[method][j,step])**2 + (Ys[method][i,step] - Ys[method][j,step])**2)
            V += MASS[i]*MASS[j]/dist
    return 0.5*T - G*V

def computeHamiltonianEvolutions():
    H = [np.zeros(nb_steps+1) for i in range(nb_methods)]
    for m in range(nb_methods):
        if useMethods[m]:
            for i in range(nb_steps+1):
                H[m][i] = computeHamiltonianAt(m,i)
    return H

Hs = computeHamiltonianEvolutions()
Hmax = max([max(Hs[m]) for m in range(nb_methods) if useMethods[m]])
Hmin = min([min(Hs[m]) for m in range(nb_methods) if useMethods[m]])

deltaAlpha = 1/trLen
T = [t0 + i*(tEnd-t0)/(nb_steps*T_year) for i in range(nb_steps+1)]

fig, (ax1, ax2) = plt.subplots(1,2)

H0_defined = False
for m in range(nb_methods):
    if useMethods[m]:
        ax1.plot(T, Hs[m], color=methodColors[m], label=methods[m])
        if not H0_defined:
            H0 = Hs[m][0]
            H0_defined = True
ax1.plot([t0,tEnd/T_year], [H0,H0], color='k', linestyle='--', label=r'$H_0$')
ax1.legend()
ax1.set_xlabel("Time (years)")
ax1.set_ylabel("Hamiltonian")

for t in range(nb_steps+1):
    plt.show(block=False)
    line = ax1.plot([T[t],T[t]], [Hmax,Hmin], color='k', linestyle='--')

    b = t+1
    a = b - trLen

    if a < 0:
        a = 0

    ax2.cla()
    ax2.set_xlabel("x (normalized)")
    ax2.set_ylabel("y (normalized)")

    for m in range(nb_methods):
        if useMethods[m]:
            for i in range(nb_planets):
                if i == 0:
                    ax2.scatter(Xs[m][i,t]/AU, Ys[m][i,t]/AU, color=methodColors[m], label=methods[m])
                else:
                    ax2.scatter(Xs[m][i,t]/AU, Ys[m][i,t]/AU, color=methodColors[m])
                X_plot = [x/AU for x in Xs[m][i,a:b]]
                Y_plot = [y/AU for y in Ys[m][i,a:b]]
                for j in range(b-a):
                    ax2.plot(X_plot[j:j+2], Y_plot[j:j+2],color=methodColors[m], alpha=(j+1)*deltaAlpha)

    ax2.legend()
    ax2.grid()
    plt.draw()
    plt.pause(0.5)
    line.pop(0).remove()
plt.show()

print("--------------------------------------------------------END-")
