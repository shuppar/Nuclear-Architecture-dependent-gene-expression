#!/usr/bin/python


import numpy as np
import sys
import math
import csv
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mp
from matplotlib.patches import Circle, Wedge, Polygon
np.set_printoptions(threshold=sys.maxsize)


##########################
# loading 48 post control files
A1 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/48hr_cont_postEdU/Data/test/G1.txt')
A2 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/48hr_cont_postEdU/Data/test/S.txt')
A3 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/48hr_cont_postEdU/Data/test/G2_M.txt')

##########################
# loading files for next high dose 48_pre J/m2

B1 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/24hr_UV_preEdU/Data/test/G1.txt')
B2 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/24hr_UV_preEdU/Data/test/S.txt')
B3 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/24hr_UV_preEdU/Data/test/G2_M.txt')

##########################
# loading files for next high dose 24_pre J/m2

C1 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/48hr_UV_preEdU/Data/test/G1.txt')
C2 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/48hr_UV_preEdU/Data/test/S.txt')
C3 = np.loadtxt('/media/shuppar/FINYEAR/2019/EU_EdU_HPG_Puro/Sphase_CellCycle_DDR/August17_A549_H2AX_53BP1_Post_EdU_48Hr_recovery/48hr_UV_preEdU/Data/test/G2_M.txt')

###################################################################################

A11 = (A1[:, 19]); A21 = (A2[:, 19]); A31 = (A3[:, 19]);
l11 = len(A11); l21 = len(A21); l31 = len(A31);
A11 = pd.DataFrame({"DSBs":A11, "Phase":l11*[1], "Which":l11*[1]});
A21 = pd.DataFrame({"DSBs":A21, "Phase":l21*[2], "Which":l21*[1]});
A31 = pd.DataFrame({"DSBs":A31, "Phase":l31*[3], "Which":l31*[1]});
A = np.concatenate((A11, A21, A31), axis=0);
# print(A)

B11 = (B1[:, 19]); B21 = (B2[:, 19]); B31 = (B3[:, 19]);
l11 = len(B11); l21 = len(B21); l31 = len(B31);
B11 = pd.DataFrame({"DSBs":B11, "Phase":l11*[1], "Which":l11*[2]});
B21 = pd.DataFrame({"DSBs":B21, "Phase":l21*[2], "Which":l21*[2]});
B31 = pd.DataFrame({"DSBs":B31, "Phase":l31*[3], "Which":l31*[2]});
B = np.concatenate((B11, B21, B31), axis=0);
# print(B)

C11 = (C1[:, 19]); C21 = (C2[:, 19]); C31 = (C3[:, 19]);
l11 = len(C11); l21 = len(C21); l31 = len(C31);
C11 = pd.DataFrame({"DSBs":C11, "Phase":l11*[1], "Which":l11*[3]});
C21 = pd.DataFrame({"DSBs":C21, "Phase":l21*[2], "Which":l21*[3]});
C31 = pd.DataFrame({"DSBs":C31, "Phase":l31*[3], "Which":l31*[3]});
C = np.concatenate((C11, C21, C31), axis=0);
# print(C)

dp = np.concatenate((A, B, C), axis=0);
dp = pd.DataFrame({"DSBs":dp[:,0], "Phase":dp[:,1], "Which":dp[:,2]});
# print(dp.Phase)
x = dp.DSBs; y = dp.Phase; z = dp.Which
print(np.min(dp.DSBs))

# fig2 = plt.figure()
# plt.plot(x, y)

# For Box Pots
ax = sns.catplot(x="Phase", y="DSBs", hue="Which",
            kind="box", data=dp);

# for Violin plots

ax = sns.catplot(x="Phase", y="DSBs", hue="Which",
            kind="violin", data=dp);            
         
# ax.set(ylim=(0, 20))
plt.show()



#########################################################################
