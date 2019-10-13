#!/usr/bin/python


import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import matplotlib.patches as mp
import csv


A, X, Y, Z = [], [], [], []

  
A = np.loadtxt('finStat.txt')
A = A[A[:,0].argsort()]				# Sorting A wrt increasing DNA Content (column 0)
X = np.loadtxt('G1.txt')
Y = np.loadtxt('S.txt')
Z = np.loadtxt('G2_M.txt')


Bf = len(X[:,0])/len(A[:,0])
# print(Bf)
Cf = len(Y[:,0])/len(A[:,0])
Df = len(Z[:,0])/len(A[:,0])
Bf = 2*np.pi*Bf
Cf = 2*np.pi*Cf
Df = 2*np.pi*Df

B = np.concatenate((X, Y, Z), axis=0)

X1 = 1 - B[:, 17]			# Col No.: for 3d: 17,18,19,20;  for 2d: 32, 33, 34, 35 
# Y1 = B[:, 38]/np.average(B[:,38])		# Protein column: 38 (size of spots)
Y1 = B[:, 36]		# mRNA column: 36 in case (size of spots)
Z1 = B[:, 21]				# mRNA column (36) or expression (0? 1?) column (21, 22, 23)  (Color of spots)

r1 = 1
t1 = 0 
t2 = Bf
t3 = Bf+Cf

# print(t1); print(t2); print(t3);

N = len(B[:,1])
# print(N)
# r = 2 * np.random.rand(N)
# theta = 2 * np.pi * np.random.rand(N)
r = X1
theta = np.linspace(start=0, stop = 2*np.pi, num = N)
area = 2*Y1		# for protein, mRNA, cell cycle plots
# area = 50 			# for allele expression, binary plots
colors = Z1

# Modified colourmap

top = [255/255, 203/255, 5/255]
bottom = [237/255, 0/255, 3/255]
newcolors = np.vstack((top, bottom))
newcmp = ListedColormap(newcolors, name='OrangeBlue')

#


fig1 = plt.figure()
ax1 = fig1.add_subplot(111, projection='polar')
cb = ax1.scatter(theta, r, c=colors, s=area, cmap=newcmp, alpha=0.75)

ax1.plot([t1, t1], [0, r1], c = [254/255,246/255,0/255], linewidth = 3, alpha = 1)
ax1.plot([t2, t2], [0, r1], c = [250/255,155/255,0/255], linewidth = 3,  alpha = 1)
ax1.plot([t3, t3], [0, r1], c = [150/255,40/255,20/255], linewidth = 3,  alpha = 1)


ax1.set_rmax(1)
ax1.set_rticks([0.8])  # Less radial ticks
plt.tick_params(labelsize=0.5, pad=6, labelrotation = 60);
ax1.set_xticks([])  # Less radial ticks
ax1.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax1.grid(True)


# cbar = fig1.colorbar(cb, ax=ax1, ticks = [0, 1])
# cb.set_clim(0,1)
plt.grid(b=None, which='major', color = [80/255,17/255,6/255], linewidth = 3, alpha = 0.5)


# ax.set_thetamin(45)
# ax.set_thetamax(135)
plt.show()
