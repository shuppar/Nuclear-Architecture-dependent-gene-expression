#!/usr/bin/python


import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.patches as mp
import csv


A, X, Y, Z = [], [], [], []

A = np.loadtxt('Intensity.txt')
X = np.loadtxt('G1.txt')
Y = np.loadtxt('S.txt')
Z = np.loadtxt('G2_M.txt')
g2 = np.loadtxt('G1_peak.txt')
g2 = g2[2]
X = X[X[:,0].argsort()]				# Sorting X wrt increasing DNA Content (column 0)
Y = Y[Y[:,0].argsort()]
Z = Z[Z[:,0].argsort()]

############################################################################################
# Normalizing the things

B = np.concatenate((X, Y, Z), axis=0)
B[:,0] = (2*np.pi*B[:,0])/g2
X[:,0] = (2*np.pi*X[:,0])/g2
Y[:,0] = (2*np.pi*Y[:,0])/g2
Z[:,0] = (2*np.pi*Z[:,0])/g2

##############################################################################################
# Getting fractions and setting plot parameters

Bf = len(X[:,0])/len(A[:,0])
Cf = len(Y[:,0])/len(A[:,0])
Df = len(Z[:,0])/len(A[:,0])
Bf = 2*np.pi*Bf
Cf = 2*np.pi*Cf
Df = 2*np.pi*Df
# r1 = np.max(B[:,19])*0.7
r1 = 20
roff = -15
t1 = 0 
t2 = Bf
t3 = Bf+Cf

##########################################################################################
# Preparing data for plots

X1 = X[:, 19]
Y1 = X[:, 2]/np.average(X[:,2])
area1 = 50*Y1
Z1 = X[:, 2]
theta1 = np.linspace(start=t1, stop = t2, num = len(X1))
# theta1 = X[:,0]

X2 = Y[:, 19]
Y2 = Y[:, 2]/np.average(Y[:,2])
area2 = 50*Y2
Z2 = Y[:, 2]
theta2 = np.linspace(start=t2, stop = t3, num = len(X2))
# theta2 = Y[:,0]

X3 = Z[:, 19]
Y3 = Z[:, 2]/np.average(Z[:,2])
area3 = 50*Y2
Z3 = Z[:, 2]
theta3 = np.linspace(start=t3, stop = 2*np.pi, num = len(X3))
# theta3 = Z[:,0]

########################################################################################
# Plotting scatters
fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='polar')
ax1.scatter(theta1, X1, c=[254/255,246/255,0/255], s=area1, edgecolors=[254/255,216/255,0/255], alpha=0.85)
ax1.scatter(theta2, X2, c=[250/255,155/255,0/255], s=area1, edgecolors=[235/255,110/255,0/255], alpha=0.85)
ax1.scatter(theta3, X3, c=[150/255,40/255,20/255], s=area1, edgecolors=[80/255,17/255,6/255], alpha=0.85)

# Positions of R0 and theta0

ax1.set_rorigin(roff)
# ax1.set_theta_zero_location('W', offset=0)

#############################################################################################
# Plotting mean level arcs
res = 1000
m1 = np.average(X1); m2 = np.average(X2); m3 = np.average(X3)
m1a = m1*np.ones(res); m2a = m2*np.ones(res); m3a = m3*np.ones(res);

e1 = np.std(X1)/np.sqrt(len(X1)); e2 = np.std(X2)/np.sqrt(len(X2)); e3 = np.std(X3)/np.sqrt(len(X3));
e1 = [m1, e1]; e2 = [m2, e2]; e3 = [m3, e3];

mt1 = np.linspace(start=t1, stop = t2, num = res)
mt2 = np.linspace(start=t2, stop = t3, num = res)
mt3 = np.linspace(start=t3, stop = 2*np.pi, num = res)

# Arcs

ax1.plot(mt1, m1a, c = [0, 0, 0], linewidth = 3.5, alpha = 1)
ax1.plot(mt2, m2a, c = [0, 0, 0], linewidth = 3.5,  alpha = 1)
ax1.plot(mt3, m3a, c = [0, 0, 0], linewidth = 3.5,  alpha = 1)

# Errors

ax1.plot([0, 0], e1 , c = [0, 0, 0], linewidth = 3.5, alpha = 1)
ax1.plot([t2, t2], e2, c = [0, 0, 0], linewidth = 3.5,  alpha = 1)
ax1.plot([t3, t3], e3, c = [0, 0, 0], linewidth = 3.5,  alpha = 1)

######################################################################################
# Plotting cell cycle phase lines

ax1.plot([t1, t1], [roff, r1], c = [254/255,246/255,0/255], linewidth = 5, alpha = 1)
ax1.plot([t2, t2], [roff, r1], c = [250/255,155/255,0/255], linewidth = 5,  alpha = 1)
ax1.plot([t3, t3], [roff, r1], c = [150/255,40/255,20/255], linewidth = 5,  alpha = 1)

#################################################################################################
# Axes Properties (always to be set in the end)

ax1.set_rmax(r1)
ax1.set_rticks([0, r1])  # Less radial ticks
ax1.set_xticks([])  # Less radial ticks
ax1.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax1.grid(True)
plt.tick_params(labelsize=0.5, pad=6, labelrotation = 60);
plt.grid(b=None, which='major', color = [80/255,17/255,6/255], linewidth = 3, alpha = 0.7)


# ax.set_thetamin(45)
# ax.set_thetamax(135)
plt.show()
