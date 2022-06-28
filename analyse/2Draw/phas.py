# -*- coding: utf-8 -*-
"""
Created on Fri May  6 08:01:23 2022

@author: 姜高晓
"""

import numpy as np
import matplotlib.pyplot as plt

W=[]

for i in range(58,20,-2):
    try:
        temp=np.loadtxt("W%d.txt"%(i))
        #W.append(temp/temp[0])
        W.append(temp)
    except:
        pass
plt.rcParams.update({"font.size":22})
plt.imshow(W,interpolation=("gaussian"),vmin=0.001)
#plt.imshow(W,interpolation=("none"),cmap="rainbow")
#W=np.random.normal((30,30),20)
plt.imshow(W,interpolation=("gaussian"),cmap="rainbow")
#plt.imshow(W,interpolation=("none"),cmap="rainbow")
plt.xlabel("$\lambda$")
plt.ylabel("$v_{active}$")
plt.xticks(np.arange(0,len(W[0]),5),np.arange(0,len(W[0]),5)*0.15)
plt.yticks(np.arange(0,len(W),4),np.arange(len(W)+10,0,-4)*2)
#plt.colorbar(label="$\omega$")
#plt.colorbar(label="$\omega/\omega_{max}$")
plt.colorbar(label="$\omega/rad\cdot τ ^{-1}$")

plt.tight_layout()
plt.savefig("xiangtu_w.svg")
#plt.savefig("xiangtu_w_w.svg")
plt.show()