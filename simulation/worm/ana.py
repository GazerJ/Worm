
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 19:08:05 2021

@author: 94917
"""

import numpy as np
import matplotlib.pyplot as plt
class read:
    def __init__(self,i,n,m):
        self.n=n
        self.dt=1
        self.strn="{:.8f}".format(n)
        self.strm="{:.8f}".format(m)
        self.stri="{:d}".format(i)
        print(self.strn,self.strm)
        self.x1=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"x1.txt")
        self.x2=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"x2.txt")
        self.x3=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"x3.txt")
        self.x4=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"x4.txt")
        self.y1=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"y1.txt")
        self.y2=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"y2.txt")
        self.y3=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"y3.txt")
        self.y4=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"y4.txt")
    def getMSD(self):
        self.lenge=[]    
        for i in range(1,len(self.x1)-1):
            self.lenge.append(((self.x1[i,:]**2+self.y1[i,:]**2+self.x2[i,:]**2+self.y2[i,:]**2+self.x3[i,:]**2+self.y3[i,:]**2+self.x4[i,:]**2+self.y4[i,:]**2).mean()-(self.x1[i,:].mean()**2+self.y1[i,:].mean()**2+self.x2[i,:].mean()**2+self.y2[i,:].mean()**2+self.x3[i,:].mean()**2+self.y3[i,:].mean()**2+self.x4[i,:].mean()**2+self.y4[i,:].mean()**2))/4)#集群的转动惯量
        plt.plot(np.arange(1,len(self.lenge)+1)*self.dt,self.lenge,label="T:"+self.strn)
        plt.yscale('log')
        plt.xscale('log')
        plt.legend(loc="upper left")
        plt.xlabel("$t/τ$")
        plt.ylabel("$MSD$")
        plt.title("$MSD$")
        self.lowMSD=np.min(self.lenge)
    def getDot(self):
        dot=[]
        for i in range(0,len(self.x1),1):
            P=np.array([self.x1[i,:]-self.x4[i,:],self.y1[i,:]-self.y4[i,:]])
            R=np.array([(self.x4[i,:]+self.x1[i,:])/2-self.x1[i,:].mean(),(self.y4[i,:]+self.y1[i,:])/2-self.y1[i,:].mean()])
            Pmod=np.sqrt(P[0,:]**2+P[1,:]**2)
            Rmod=np.sqrt(R[0,:]**2+R[1,:]**2)
            P=P/Pmod
            R=R/Rmod
            dot.append((((P*R)[0,:]+(P*R)[1,:])).mean())
        plt.plot(np.arange(len(dot))*self.dt*100,dot,label="T:"+self.strn)
        plt.title("$P·R$")
        plt.xlabel("$t/τ$")
        plt.ylim(-1, 1)
        plt.legend(loc="upper left")
        self.dotmean=np.array(dot)[-200:].mean()
    def getCross(self):
        cross=[]
        for i in range(0,len(self.x1),1):
            P=np.array([self.x1[i,:]-self.x4[i,:],self.y1[i,:]-self.y4[i,:]])
            R=np.array([(self.x4[i,:]+self.x1[i,:])/2-self.x1[i,:].mean(),(self.y4[i,:]+self.y1[i,:])/2-self.y1[i,:].mean()])
            Pmod=np.sqrt(P[0,:]**2+P[1,:]**2)
            Rmod=np.sqrt(R[0,:]**2+R[1,:]**2)
            P=P/Pmod
            R=R/Rmod
            per=[]
            for j in range(len(P)):
                per.append(np.abs(np.linalg.det([P[:,j],R[:,j]])))
            cross.append((np.array(per)).mean())
        plt.plot(np.arange(len(cross))*self.dt*100,cross,label="T:"+self.strm)
        plt.ylim(-1, 1)
        plt.xlabel("$t/τ$")
        plt.title("$PXR$")
        plt.legend(loc="upper left")
        self.crossmean=np.array(cross)[-200:].mean()
    def getFig(self):
        for i in range(0,len(self.x1),10):
    #        for i in range(0,10000,1000):
            for j in range(len(self.x1[0,:])):
                plt.plot([self.x1[i,j],self.x2[i,j],self.x3[i,j],self.x4[i,j]],[self.y1[i,j],self.y2[i,j],self.y3[i,j],self.y4[i,j]],lw=4)
                plt.scatter(self.x1[i,j],self.y1[i,j])
            plt.xlim(-30,30)
            plt.ylim(-30,30)
            plt.show()
    def getCrow(self):
        MSD=(np.log10(self.lenge))
        crow=np.where(MSD<2*MSD.min())[0]
        tau=(crow.max()-crow.min())*self.dt
        return MSD.min(),tau
box=[]

m=np.arange(16)*0.00000001
for i in range(1):
    for j in range(len(m)):
        box.append(read(i,m[j],0))
        box[-1].getMSD()
plt.show()



cros=[]
dot=[]
lowMSD=[]
for i in box:
    i.getCross()
    cros.append(i.crossmean)
plt.ylim(-1,1)
plt.show()
for i in box:
    i.getDot()
    dot.append(i.dotmean)
    lowMSD.append(i.lowMSD)
plt.ylim(-1,1)
plt.show()
'''
m=list(m)
mm=m+m+m
plt.scatter(mm,dot)
plt.ylim(-1,1)
plt.show()

plt.scatter(mm,cros)
plt.ylim(-1,1)
plt.show()

plt.scatter(mm,lowMSD)
plt.show()

n=20
dotMean=np.zeros(n)
crosMean=np.zeros(n)
lowMSDMean=np.zeros(n)
for i in range(n):
    dotMean[i]=(box[i].dotmean+box[i+n].dotmean+box[i+2*n].dotmean)/3
    crosMean[i]=(box[i].crossmean+box[i+n].crossmean+box[i+2*n].crossmean)/3
    lowMSDMean[i]=(box[i].lowMSD+box[i+n].lowMSD+box[i+2*n].lowMSD)/3

plt.plot(m,dotMean)
plt.ylim(-1,1)
plt.show()
plt.plot(m,crosMean)
plt.ylim(-1,1)
plt.show()
plt.plot(m,lowMSDMean)
plt.show()


'''
'''
t=[]
Tsun=[]
for n in T:
    for m in Ts:
        Tsun.append(n)    
for i in box:
    x=i.getCrow()
    t.append(x[1])
plt.plot(Tsun,t)
plt.yscale('log')

plt.show()

'''
