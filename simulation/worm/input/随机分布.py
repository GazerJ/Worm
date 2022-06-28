# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 09:02:43 2021

@author: 94917
"""
import matplotlib.pyplot as plt
import numpy as np 
r0=1
N=1000
L=1000
x4=np.random.random(N)*L
y4=np.random.random(N)*L
s=np.random.random(N)*6.28
x3=x4+r0*np.cos(s)
y3=y4+r0*np.sin(s)

x2=x4+2*r0*np.cos(s)
y2=y4+2*r0*np.sin(s)

x1=x4+3*r0*np.cos(s)
y1=y4+3*r0*np.sin(s)


def x_y_(x,y,s):
    return np.cos(s)*x-np.sin(s)*y,np.sin(s)*x+np.cos(s)*y
#第一层

for j in range(len(x1)):
        plt.plot([x1[j],x2[j],x3[j],x4[j]],[y1[j],y2[j],y3[j],y4[j]])
        plt.scatter(x1[j],y1[j])
plt.show()



np.savetxt("toheart30x1.txt",x1)
np.savetxt("toheart30x2.txt",x2)
np.savetxt("toheart30x3.txt",x3)
np.savetxt("toheart30x4.txt",x4)
np.savetxt("toheart30y1.txt",y1)
np.savetxt("toheart30y2.txt",y2)
np.savetxt("toheart30y3.txt",y3)
np.savetxt("toheart30y4.txt",y4)