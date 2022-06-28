# -*- coding: utf-8 -*-
"""
Created on Sat Nov 20 10:07:07 2021

@author: GazerJ
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Nov 14 09:02:43 2021

@author: 94917
"""
import matplotlib.pyplot as plt
import numpy as np 



def x_y_new(x,y,s):
    return np.cos(s)*x-np.sin(s)*y,np.sin(s)*x+np.cos(s)*y



# n度对称性
n=8
r0=1


N=n*6
x1_std=N-13
x2_std=r0+x1_std
x3_std=r0+x2_std
x4_std=r0+x3_std
y1_std=0
y2_std=0
y3_std=0
y4_std=0
s=np.arange(0,2*np.pi,2*np.pi/N)
x1new0,y1new0=x_y_new(x1_std,y1_std,s)
x2new0,y2new0=x_y_new(x2_std,y2_std,s)
x3new0,y3new0=x_y_new(x3_std,y3_std,s)
x4new0,y4new0=x_y_new(x4_std,y4_std,s)



N=n*5
x1_std=N-12
x2_std=r0+x1_std
x3_std=r0+x2_std
x4_std=r0+x3_std
y1_std=0
y2_std=0
y3_std=0
y4_std=0
s=np.arange(0,2*np.pi,2*np.pi/N)
x1new1,y1new1=x_y_new(x1_std,y1_std,s)
x2new1,y2new1=x_y_new(x2_std,y2_std,s)
x3new1,y3new1=x_y_new(x3_std,y3_std,s)
x4new1,y4new1=x_y_new(x4_std,y4_std,s)
'''
x1new1=np.append(x1new1,x1new0)
x2new1=np.append(x2new1,x2new0)
x3new1=np.append(x3new1,x3new0)
x4new1=np.append(x4new1,x4new0)
y1new1=np.append(y1new1,y1new0)
y2new1=np.append(y2new1,y2new0)
y3new1=np.append(y3new1,y3new0)
y4new1=np.append(y4new1,y4new0)
'''

'''
x1new1,y1new1=[],[]
x2new1,y2new1=[],[]
x3new1,y3new1=[],[]
x4new1,y4new1=[],[]
'''



N=n*4
x1_std=N-11
x2_std=r0+x1_std
x3_std=r0+x2_std
x4_std=r0+x3_std
y1_std=0
y2_std=0
y3_std=0
y4_std=0
s=np.arange(0,2*np.pi,2*np.pi/N)
x1new2,y1new2=x_y_new(x1_std,y1_std,s)
x2new2,y2new2=x_y_new(x2_std,y2_std,s)
x3new2,y3new2=x_y_new(x3_std,y3_std,s)
x4new2,y4new2=x_y_new(x4_std,y4_std,s)
'''
x1new2,y1new2=[],[]
x2new2,y2new2=[],[]
x3new2,y3new2=[],[]
x4new2,y4new2=[],[]
'''
x1new=np.append(x1new1,x1new2)
x2new=np.append(x2new1,x2new2)
x3new=np.append(x3new1,x3new2)
x4new=np.append(x4new1,x4new2)
y1new=np.append(y1new1,y1new2)
y2new=np.append(y2new1,y2new2)
y3new=np.append(y3new1,y3new2)
y4new=np.append(y4new1,y4new2)


N=n*3
x1_std=N-10
x2_std=r0+x1_std
x3_std=r0+x2_std
x4_std=r0+x3_std
y1_std=0
y2_std=0
y3_std=0
y4_std=0
s=np.arange(0,2*np.pi,2*np.pi/N)
x1new2,y1new2=x_y_new(x1_std,y1_std,s)
x2new2,y2new2=x_y_new(x2_std,y2_std,s)
x3new2,y3new2=x_y_new(x3_std,y3_std,s)
x4new2,y4new2=x_y_new(x4_std,y4_std,s)
'''
x1new2,y1new2=[],[]
x2new2,y2new2=[],[]
x3new2,y3new2=[],[]
x4new2,y4new2=[],[]
'''

x1new=np.append(x1new,x1new2)
x2new=np.append(x2new,x2new2)
x3new=np.append(x3new,x3new2)
x4new=np.append(x4new,x4new2)
y1new=np.append(y1new,y1new2)
y2new=np.append(y2new,y2new2)
y3new=np.append(y3new,y3new2)
y4new=np.append(y4new,y4new2)




N=n*2
x1_std=N-9
x2_std=r0+x1_std
x3_std=r0+x2_std
x4_std=r0+x3_std
y1_std=0
y2_std=0
y3_std=0
y4_std=0
s=np.arange(0,2*np.pi,2*np.pi/N)
x1new2,y1new2=x_y_new(x1_std,y1_std,s)
x2new2,y2new2=x_y_new(x2_std,y2_std,s)
x3new2,y3new2=x_y_new(x3_std,y3_std,s)
x4new2,y4new2=x_y_new(x4_std,y4_std,s)

x1new=np.append(x1new,x1new2)
x2new=np.append(x2new,x2new2)
x3new=np.append(x3new,x3new2)
x4new=np.append(x4new,x4new2)
y1new=np.append(y1new,y1new2)
y2new=np.append(y2new,y2new2)
y3new=np.append(y3new,y3new2)
y4new=np.append(y4new,y4new2)

N=n
x1_std=N-5  
x2_std=r0+x1_std
x3_std=r0+x2_std
x4_std=r0+x3_std
y1_std=0
y2_std=0
y3_std=0
y4_std=0
s=np.arange(0,2*np.pi,2*np.pi/N)
x1new2,y1new2=x_y_new(x1_std,y1_std,s)
x2new2,y2new2=x_y_new(x2_std,y2_std,s)
x3new2,y3new2=x_y_new(x3_std,y3_std,s)
x4new2,y4new2=x_y_new(x4_std,y4_std,s)

x1new=np.append(x1new,x1new2)+50
x2new=np.append(x2new,x2new2)+50
x3new=np.append(x3new,x3new2)+50
x4new=np.append(x4new,x4new2)+50
y1new=np.append(y1new,y1new2)+50
y2new=np.append(y2new,y2new2)+50
y3new=np.append(y3new,y3new2)+50
y4new=np.append(y4new,y4new2)+50

for j in range(len(x1new)):
        #plt.scatter([x1new[j],x2new[j],x3new[j],x4new[j]],[y1new[j],y2new[j],y3new[j],y4new[j]],c="greed")
        plt.scatter(x1new[j],y1new[j],c="red")
        plt.scatter(x2new[j],y2new[j],c="green")
        plt.scatter(x3new[j],y3new[j],c="green")
        plt.scatter(x4new[j],y4new[j],c="green")

plt.savefig("input.svg")
plt.show()




'''
s=[]

def cultheata(c,s):
    return np.sign(np.arcsin(s))*np.arccos(c)


for i in range(len(x1new)):
    s.append(cultheata((x1new[i]-x4new[i])/np.sqrt((x1new[i]-x4new[i])**2+(y1new[i]-y4new[i])**2),(y1new[i]-y4new[i])/np.sqrt((y1new[i]-y4new[i])**2+(y1new[i]-y4new[i])**2)))



np.savetxt("s.txt",s)


'''










np.savetxt("toheart30x1.txt",x1new)
np.savetxt("toheart30x2.txt",x2new)
np.savetxt("toheart30x3.txt",x3new)
np.savetxt("toheart30x4.txt",x4new)
np.savetxt("toheart30y1.txt",y1new)
np.savetxt("toheart30y2.txt",y2new)
np.savetxt("toheart30y3.txt",y3new)
np.savetxt("toheart30y4.txt",y4new)