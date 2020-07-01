print("------------------------------------------------------------")
print("-----------------   ELLIPSE PROBLEM CODE   -----------------")
print("-------   Alexandre DUTKA - ISAE-SUPAERO - 06/2020   -------")
print("------------------------------------------------------------")

import matplotlib.pyplot as plt


#%% PARAMETERS

#initial conditions
x0 = 1.
y0 = 0.

#ellipse parameters
a = 2.
b = 1.

#integration parameters
dt = 0.1
nb_steps = 200


#%% SCRIPT

X1 = [x0]
Y1 = [y0]

for t in range(nb_steps):
    X1.append(X1[-1] - dt*Y1[-1]/b**2)
    Y1.append(Y1[-1] + dt*X1[-1]/a**2)

X2 = [x0]
Y2 = [y0]

factorM = (1/a**2 - 1/b**2)
factorP = (1/a**2 + 1/b**2)

for t in range(nb_steps):
    X2.append(X2[-1]*(1+dt*factorM) - dt*factorP*Y2[-1])
    Y2.append(Y2[-1]/(1+dt*factorM) + X2[-1]*dt*factorP/(1+dt*factorM))

plt.figure()
plt.xlim(-2,2)
plt.ylim(-2,2)
plt.plot(X1,Y1, label="Separable case")
plt.plot(X2,Y2, label="Non-separable case")
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.show(block=False)

T = [i*dt for i in range(nb_steps+1)]
H1 = [0.5*((X1[i]/a)**2 + (Y1[i]/b)**2) for i in range(nb_steps+1)]
H2 = [0.5*(((X2[i]-Y2[i])/a)**2 + ((X2[i]+Y2[i])/b)**2) for i in range(nb_steps+1)]

H1n = [h/H1[0] for h in H1]
H2n = [h/H2[0] for h in H2]

plt.figure()
plt.ylim(0.,1.5)
plt.plot(T,H1n, label="Separable case")
plt.plot(T,H2n, label="Non-separable case")
plt.xlabel("time")
plt.ylabel(r"Normalized Hamiltonian $H^* = H \; / \; H_0$")
plt.show()

print("--------------------------------------------------------END-")
