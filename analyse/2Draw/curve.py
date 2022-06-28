# -*- coding: utf-8 -*-
"""
Created on Thu Apr 28 12:54:23 2022

@author: 姜高晓
"""

import numpy as np
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
plt.rcParams.update({"font.size":15})

def Draw(data,x,v,s):
    def func(x,a,b,c,d,e):
        return b*np.exp(-a*(x))
    def funccc(x,a,b,c,d,e):
        #return a/(c+ np.exp(-b*x**3))
        return v*( a+c*np.exp(-b*(x+d)*2))
    def funcc(x,a,b,c,d,e):
        return v*v/(a*x+b*v)*c   
    
    
    
    popt, pcov = curve_fit(funcc,x,data)
    
    yfit=[funcc(i,popt[0],popt[1],popt[2],popt[3],popt[4]) for i in x]
    plt.ylim(0,0.004)
    #plt.xlim(0.01,0.9)
    #plt.xscale("log");plt.yscale("log")
    #plt.plot(x,yfit,label="$"+"y=%d/(%d +%d x+%.5f)"%(v,popt[0],popt[1],popt[2])+"$",color="black")
    plt.plot(x,yfit,color="black")
    plt.xlabel("$\lambda$")
    
    plt.ylabel("$\omega$")
    
    
    plt.scatter(x,data,s=100,marker=s,label="$v_{active}=%d"%(v)+"$")
    plt.tight_layout()
    plt.legend(loc="upper right")
    plt.savefig("wlog.svg")

'''
data = np.loadtxt("W20.txt")
x=np.arange(len(data))*0.15
Draw(data,x,20)
'''

'''
data = np.loadtxt("W22.txt")
x=np.arange(len(data))*0.15
Draw(data,x,22)  
data = np.loadtxt("W24.txt")
x=np.arange(len(data))*0.15
Draw(data,x,24)  
data = np.loadtxt("W26.txt")
x=np.arange(len(data))*0.15
Draw(data,x,26)  
data = np.loadtxt("W28.txt")
x=np.arange(len(data))*0.15
Draw(data,x,28)  
'''

data = np.loadtxt("W30.txt")[1:-2:1]
x=np.arange(len(data))*0.15+0.15
Draw(data,x,30,"o")  

'''
data = np.loadtxt("W32.txt")
x=np.arange(len(data))*0.15
Draw(data,x,32)  

data = np.loadtxt("W34.txt")
x=np.arange(len(data))*0.15
Draw(data,x,34)  

data = np.loadtxt("W36.txt")
x=np.arange(len(data))*0.15
Draw(data,x,36)  
data = np.loadtxt("W38.txt")
x=np.arange(len(data))*0.15
Draw(data,x,38)  

'''
data = np.loadtxt("W40.txt")[1:-2:1]
x=np.arange(len(data))*0.15+0.15
Draw(data,x,40,"2")  
'''
data = np.loadtxt("W42.txt")
x=np.arange(len(data))*0.15
Draw(data,x,42)  

data = np.loadtxt("W44.txt")
x=np.arange(len(data))*0.15
Draw(data,x,44) 
'''

data = np.loadtxt("W50.txt")[1:-2:1]
x=np.arange(len(data))*0.15+0.15
Draw(data,x,50,"x") 


'''
data = np.loadtxt("W58.txt")
x=np.arange(len(data))*0.15
Draw(data,x,58) 

'''
'''
data = np.loadtxt("W40.txt")/2
x=np.arange(len(data))*0.02
Draw(data,x,40)
plt.show()
'''
