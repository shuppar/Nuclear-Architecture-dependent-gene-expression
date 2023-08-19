#!/usr/bin/python


import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
import matplotlib.patches as mp
import csv


A, X, Y, Z = [], [], [], []

  
A = np.loadtxt('finStat.txt')
A = A[A[:,0].argsort()]				# Sorting A wrt increasing DNA Content (column 0)
X = np.loadtxt('G1.txt')
Y = np.loadtxt('S.txt')
Z = np.loadtxt('G2_M.txt')


Bf = len(X[:,1])/len(A[:,1])
Cf = len(Y[:,1])/len(A[:,1])
Df = len(Z[:,1])/len(A[:,1])
Bf = 2*np.pi*Bf
Cf = 2*np.pi*Cf
Df = 2*np.pi*Df


B = np.concatenate((X, Y, Z), axis=0)

X1 = 1 - A[:, 7]/A[:,2]			# Changed it back to A from B after sorting A wrt increasing DNA Content.
Y1 = A[:, 14]/np.average(B[:,14])
Z1 = A[:, 12]

r1 = 1
t1 = 0 
t2 = Bf
t3 = Bf+Cf


N = len(A[:,1])
# print(N)
# r = 2 * np.random.rand(N)
# theta = 2 * np.pi * np.random.rand(N)
r = X1
theta = np.linspace(start=0, stop = 2*np.pi, num = N)
area = 50*Y1
colors = Z1


fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='polar')
c = ax1.scatter(theta, r, c=colors, s=area, cmap='plasma', alpha=0.75)
ax1.plot([0, 0], [0, r1], c = [254/255,246/255,0/255], linewidth = 3, alpha = 0.75)
ax1.plot([t2, t2], [0, r1], c = [250/255,155/255,0/255], linewidth = 3,  alpha = 0.75)
ax1.plot([t3, t3], [0, r1], c = [150/255,40/255,20/255], linewidth = 3,  alpha = 0.75)


ax1.set_rmax(1)
ax1.set_rticks([0.8])  # Less radial ticks
plt.tick_params(labelsize=0.5, pad=6, labelrotation = 60);
ax1.set_xticks([])  # Less radial ticks
ax1.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax1.grid(True)


fig1.colorbar(c, ax=ax1)
plt.grid(b=None, which='major', color = [80/255,17/255,6/255], linewidth = 3, alpha = 0.5)


# ax.set_thetamin(45)
# ax.set_thetamax(135)
plt.show()
