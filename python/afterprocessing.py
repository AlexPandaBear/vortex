print("------------------------------------------------------------")
print("------------------   VORTEX METHOD CODE   ------------------")
print("-------------------   After-Processing   -------------------")
print("-------   Alexandre DUTKA - ISAE-SUPAERO - 05/2020   -------")
print("------------------------------------------------------------")


#%% IMPORTS

print("Importing libraries")

import matplotlib.pyplot as plt
import time
import math
import numpy as np
import _vortex as vtx


#%% INSTRUCTIONS & PLOTTING PARAMETERS

print("Reading instructions")

loadNewData = True #if False, previously loaded data will be used (works only if python is used in interactive mode (-i) or in an IDE)
dataFolder = "../data"
dataFile = "test_sv"

plotVtxAnimation = False
vtxMvt_reframe = True
traceLength = 5

plotVtxConfig = False
vtxConfig_reframe = True
steps_vtxConfig = [35]#[i for i in range(101) if i%10==0]

plotCompositionField = False
steps_compositionField = [i for i in range(101) if i%10==0]
compositionRadius = "auto"
plotCompoAnimation = False

plotVelocityField = True
steps_velocityField = [55] #[i for i in range(101) if i%10==0]
showVelocityVectors = True

plotVorticityField = False
steps_vorticityField = [i for i in range(101) if i%10==0]
derivationDistance = "auto"

plotPressureCoefField = False
steps_pressureCoefField = [0, 100, 200, 500, 1000]

plotStreamlines = False
steps_streamlines = [i for i in range(101) if i%10==0]
integrationStep = 0.1
nb_streamlines = 20

plotHamiltonianEvolution = False

numberOfThreads = 4


#%% FUNCTIONS

print("Loading functions")

#Loads the necessary workspace, saved in the file file.wsdata (associated to the simulation file.simdata)
def loadWorkspace(file):
    f = open(file, "r")
    global width
    width = float(f.readline())
    global height
    height = float(f.readline())
    global U0
    U0 = float(f.readline())
    global RSLT
    RSLT = float(f.readline())
    global RIHPA
    RIHPA = float(f.readline())
    global nx
    nx = int(f.readline())
    global ny
    ny = int(f.readline())
    global t0
    t0 = float(f.readline())
    global tEnd
    tEnd = float(f.readline())
    global nb_steps
    nb_steps = int(f.readline())
    global temporalIntegrationMethod
    temporalIntegrationMethod = f.readline().rstrip()
    global periodicity
    periodicity = f.readline().rstrip()
    if periodicity == "True":
        periodicity = True
    else:
        periodicity = False
    global dx
    dx = width/(nx-1)
    global dy
    dy = height/(ny-1)
    global SLT
    SLT = 0.5*height*RSLT
    global regRadius
    regRadius = float(f.readline())
    f.close()

#Plots an animation of the fluid movement computed in the sim. manager s, with a trace length trLen for the vortices' trajectories in number of steps (min. 2)
def vtxMvtAnimation(s, trLen):
	X = [s.getVtxXs(i) for i in range(nx*ny)]
	Y = [s.getVtxYs(i) for i in range(nx*ny)]

	if vtxMvt_reframe:
		X_periodic = []
		for i in range(nx*ny):
			X_periodic.append([])
			for t in range(nb_steps):
				x = X[i][t]
				if x < 0.:
					x += width
				elif x > width:
					x -= width
				X_periodic[i].append(x)

	plt.figure()
	plt.show(block=False)

	for t in range(nb_steps):
		plt.cla()

		b = t+1
		a = b - trLen

		if a < 0:
			a = 0

		for i in range(nx*ny):
			colorSwitch = (Y[i][0] > 0.)
			color = 'b'

			if colorSwitch:
				color = 'r'

			if vtxMvt_reframe:
				X_plot = X_periodic[i][a:b]
			else:
				X_plot = X[i][a:b]

			Y_plot = Y[i][a:b]

			if vtxMvt_reframe:
				plt.xlim(0, width)
				plt.ylim(-0.6*height, 0.6*height)
				if (max(X_plot)-min(X_plot) < 0.9*width):
					plt.plot(X_plot, Y_plot, color)
			else:
				plt.plot(X_plot, Y_plot, color)

		plt.draw()
		plt.pause(1e-17)
		time.sleep(0.05)
	plt.show()

#Plots the positions of the vortices (with respective circulations represented in the dots' sizes), along with the velocity and vorticity profiles at step step
def vtxConfig(s, step):
    X = s.getXsAt(step)
    Y = s.getYsAt(step)
    C = s.getCirculationsAt(step)

    nb_vtx = len(X)

    if vtxConfig_reframe:
    	for i in range(nb_vtx):
    		if X[i] > width:
    			X[i] -= width
    		elif X[i] < 0.:
    			X[i] += width

    circMax = max([abs(circ) for circ in C])

    plt.figure()
    plt.title("Vortex configuration at step {} (ie t = {} sec)".format(step, T[step]))

    for i in range(nb_vtx):
    	color = 'dodgerblue'
    	if C[i] > 0.:
    		color = 'hotpink'
    	plt.scatter(X[i], Y[i], 100*abs(C[i])/circMax, c=color)

    plt.xlabel("x (m)")
    plt.ylabel("y (m)")
    plt.ylim([-0.5*height, 0.5*height])

    plt.show()

#Plots the fluid composition at the step step in the sim. manager s, using a matrix of nbx*nby points and considering vortices inside a radius of radius to compute
def fluidComposition(s, nbx, nby, step, radius):
    plt.figure()
    plt.title("Composition field at step {} (ie t = {} sec)".format(step, T[step]))

    x = np.linspace(0, width, nbx)
    y = np.linspace(-0.5*height, 0.5*height, nby)
    X,Y = np.meshgrid(x,y)
    C = np.zeros((nby,nbx))

    for i in range(nbx):
    	for j in range(nby):
            compo = s.computeCompositionAt(x[i], y[j], step, radius)
            C[j,i] = compo[0] - compo[1]
            if math.isnan(C[j,i]):
            	C[j,i] = 0.

    plt.imshow(C, extent=[0, width, -0.5*height, 0.5*height], origin="lower", alpha=.75, cmap='seismic')
    #c = plt.contour(X, Y, C, colors='black')
    #plt.clabel(c)

    plt.show()

#Plots an animation of the evolution of the composition field of the simulation of sim. manager s, using a matrix of nbx*nby points and considering vortices inside a radius of radius to compute
def fluidCompoAnim(s, nbx, nby, radius):
    x = np.linspace(0, width, nbx)
    y = np.linspace(-0.5*height, 0.5*height, nby)
    X,Y = np.meshgrid(x,y)
    C = np.zeros((nby,nbx))

    plt.figure()
    plt.show(block=False)

    for step in range(nb_steps+1):
        plt.cla()

        for i in range(nbx):
            for j in range(nby):
                compo = s.computeCompositionAt(x[i], y[j], step, radius)
                C[j,i] = compo[0] - compo[1]
                if math.isnan(C[j,i]):
                	C[j,i] = 0.

        plt.contourf(X, Y, C, alpha=.75, cmap='seismic')
        c = plt.contour(X, Y, C, colors='black')
        plt.clabel(c)

        plt.annotate("Step {}/{}".format(step, nb_steps), xy=(0.02*width,0.45*height))

        plt.draw()
        plt.pause(1e-17)
        time.sleep(0.05)

    plt.show()

#Plots the velocity field at the step step in the sim. manager s, using a matrix of nbx*nby points to compute
def velocityField(s, nbx, nby, step):
    plt.figure()
    plt.title("Velocity field at step {} (ie t = {} sec)".format(step, T[step]))
    x = np.linspace(0, width, nbx)
    y = np.linspace(-0.5*height, 0.5*height, nby)

    X,Y = np.meshgrid(x,y)

    U = np.zeros((nby,nbx))
    V = np.zeros((nby,nbx))
    M = np.zeros((nby,nbx))
    A = np.zeros((nby,nbx))

    for i in range(nbx):
    	for j in range(nby):
            res = s.computeVelocityAt(x[i], y[j], step, periodicity, width)
            U[j,i] = res[0]
            V[j,i] = res[1]
            M[j,i] = res[2]
            A[j,i] = res[3]

    h = min(width/(nbx-1), height/(nby-1))

    if showVelocityVectors:
        for i in range(nbx):
        	for j in range(nby):
        		mag_star = M[j,i]/U0
        		arg = A[j,i]
        		plt.annotate('', xy=(x[i]+h*mag_star*np.cos(arg),y[j]+h*mag_star*np.sin(arg)), xytext=(x[i],y[j]), arrowprops={'arrowstyle': '->', 'lw': 2, 'color': 'white'}, va='center')

    plt.contourf(X, Y, M, alpha=.75, cmap='jet')
    c = plt.contour(X, Y, M, colors='black')
    plt.clabel(c)

    plt.show()

#Plots the vorticity field's maggnitude at step step computed in sim. manager s, using nbx points of calculus horizontaly, and nby verticaly. Numerical derivation is done with staptial step h
def vorticityContours(s, nbx, nby, step, h):
    if h == "auto":
        h = max(dx, dy)

    plt.figure()
    plt.title("Vorticity field at step {} (ie t = {} sec)".format(step, T[step]))
    x = np.linspace(0, width, nbx)
    y = np.linspace(-0.5*height, 0.5*height, nby)

    X,Y = np.meshgrid(x,y)

    O = np.zeros((nby, nbx))

    for i in range(nbx):
    	for j in range(nby):
    		O[j,i] = sm.computeVorticityAt(x[i], y[j], step, h, periodicity, width)

    plt.contourf(X, Y, O, alpha=.75, cmap='jet')
    c = plt.contour(X, Y, O, colors='black')
    plt.clabel(c)

    plt.show()

#Plots the streamlines at step step of the simulation s, using nb streamlines and integrating each one with a h spatial step
def streamlines(s, nb, step, h):
    plt.figure()
    plt.title("Streamlines at step {} (ie t = {} sec)".format(step, T[step]))
    Y0 = np.linspace(-0.5*height, 0.5*height, nb)
    pts = [(0,y) for y in Y0 if y*U0 > 0] + [(width,y) for y in Y0 if y*U0 < 0]
    vtx_height = height

    for point in pts:
        X = [point[0]]
        Y = [point[1]]
        while (0 <= X[-1] <= width) and (-0.5*height <= Y[-1] <= 0.5*height):
            res = s.computeVelocityAt(X[-1], Y[-1], step, periodicity, width)
            X.append(X[-1] + h*res[0]/res[2])
            Y.append(Y[-1] + h*res[1]/res[2])
        plt.plot(X, Y, c='royalblue')

        yMax = max([abs(y) for y in Y])
        if yMax < vtx_height:
            vtx_height = yMax

    yLim = vtx_height - (height-vtx_height)/nb
    pts = [(0.5*width,0.5*y) for y in Y0 if 0 < 0.5*y < yLim]

    for point in pts:
        X = [point[0]]
        Y = [point[1]]
        while True:
            res = s.computeVelocityAt(X[-1], Y[-1], step, periodicity, width)
            newX = X[-1] + 0.1*h*res[0]/res[2]
            newY = Y[-1] + 0.1*h*res[1]/res[2]
            if not((0 <= newX <= width)) or not((-0.5*height <= newY <= 0.5*height)):
                X = []
                Y = []
                break
            elif newY > 0 and abs(newX-X[0]) < h and len(X) > 10:
                X.append(X[0])
                Y.append(Y[0])
                break;
            X.append(newX)
            Y.append(newY)
        plt.plot(X, Y, c='royalblue')

    plt.show()

#Plots the pressure coefficient field at the step step in the sim. manager s, using a matrix of nbx*nby points to compute
def pressureCoef(s, nbx, nby, step):
    plt.figure()
    plt.title("Pressure coefficient at step {} (ie t = {} sec)".format(step, T[step]))
    x = np.linspace(0, width, nbx)
    y = np.linspace(-0.5*height, 0.5*height, nby)

    X,Y = np.meshgrid(x,y)

    U = np.zeros((nby,nbx))
    V = np.zeros((nby,nbx))
    M = np.zeros((nby,nbx))
    A = np.zeros((nby,nbx))
    P = np.zeros((nby,nbx))

    for i in range(nbx):
    	for j in range(nby):
            res = s.computeVelocityAt(x[i], y[j], step, periodicity, width)
            U[j,i] = res[0]
            V[j,i] = res[1]
            M[j,i] = res[2]
            A[j,i] = res[3]
            P[j,i] = 1-(M[j,i]/U0)**2

    plt.contourf(X, Y, P, alpha=.75, cmap='jet')
    c = plt.contour(X, Y, P, colors='black')
    plt.clabel(c)

    plt.show()

#Plots the evolution of the Hamiltonian of the system during the simulation stored in the sim. manager s
def hamiltonianEvolution(s, nbThreads):
    plt.figure()
    plt.title("Evolution of the Hamiltonian over the simulation")

    plt.subplot(2,1,1)
    H = s.computeHamiltonianEvolution(nbThreads)
    plt.plot(T, H)
    plt.ylabel("Hamiltonian")

    plt.subplot(2,1,2)
    E = [h/H[0] - 1 for h in H]
    plt.plot(T, E)
    plt.xlabel("time (sec)")
    plt.ylabel("Relative error on the Hamiltonian")

    plt.show()


#%% SRCIPT

if loadNewData:
    if dataFile == "":
        print("Unable to load file (bad simLoadFile name)")
    else:
        print("Loading simulation data from file " + dataFolder + "/" + dataFile + ".simdata")
        sm = vtx.SM()
        sm.loadSim(dataFolder + "/" + dataFile + ".simdata")
        sm.setName(dataFile)

        print("Loading workspace data from file " + dataFolder + "/" + dataFile + ".wsdata")
        loadWorkspace(dataFolder + "/" + dataFile + ".wsdata")

        print("Data loaded successfully")

T = sm.getTimeList()

if plotVtxAnimation:
	print("Plotting vortex animation...")
	vtxMvtAnimation(sm, traceLength)

if plotVtxConfig :
	for step in steps_vtxConfig:
		print("Plotting vortices configuration at step {} (ie t = {} sec)".format(step, T[step]))
		vtxConfig(sm, step)

if compositionRadius == "auto":
    compositionRadius = 3*max(dx, dy)

if plotCompositionField:
	for step in steps_compositionField:
		print("Plotting composition field at step {} (ie t = {} sec)".format(step, T[step]))
		fluidComposition(sm, nx, ny, step, compositionRadius)

if plotCompoAnimation:
    print("Plotting composition animation...")
    fluidCompoAnim(sm, nx, ny, compositionRadius)

if plotVelocityField:
	for step in steps_velocityField:
		print("Plotting velocity field at step {} (ie t = {} sec)".format(step, T[step]))
		velocityField(sm, nx, ny, step)

if plotVorticityField:
	for step in steps_vorticityField:
		print("Plotting vorticity field at step {} (ie t = {} sec)".format(step, T[step]))
		vorticityContours(sm, nx, ny, step, derivationDistance)

if plotPressureCoefField:
    for step in steps_pressureCoefField:
        print("Plotting pressure coefficient at step {} (ie t = {} sec)".format(step, T[step]))
        pressureCoef(sm, nx, ny, step)

if plotStreamlines:
    for step in steps_streamlines:
        print("Plotting streamlines at step {} (ie t = {} sec)".format(step, T[step]))
        streamlines(sm, nb_streamlines, step, integrationStep)

if plotHamiltonianEvolution:
	print("Plotting Hamiltonian evolution")
	hamiltonianEvolution(sm, numberOfThreads)

print("--------------------------------------------------------END-")
