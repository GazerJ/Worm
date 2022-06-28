
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 19:08:05 2021

@author: 94917
"""
import sys
sys.path.append('./')
#import jplt
import numpy as np
import matplotlib.pyplot as plt
import threading
from sklearn.cluster import DBSCAN
class read(threading.Thread):
    def __init__(self,i,n,m,v):
        threading.Thread.__init__(self)
        l=30
        self.n=n
        self.dt=0.1*l
        lam=3
        self.strn="{:.8f}".format(n)
        self.strm="{:.8f}".format(m)
        self.strv="{:.1f}".format(v)
        self.id=i
        self.stri="{:d}".format(i)
        print('ID'+self.stri+'T:'+self.strn+'Ts:'+self.strm+'v0:'+self.strv)
        self.x1=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"x1.txt")[::l]
        self.x2=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"x2.txt")[::l]
        self.x3=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"x3.txt")[::l]
        self.x4=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"x4.txt")[::l]
        self.y1=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"y1.txt")[::l]
        self.y2=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"y2.txt")[::l]
        self.y3=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"y3.txt")[::l]
        self.y4=np.loadtxt("./output/data"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+"y4.txt")[::l]

        self.vx=np.zeros((len(self.x1)-1,len(self.x1[0,:])))
        self.vy=np.zeros((len(self.y1)-1,len(self.y1[0,:])))
        self.xmean=np.zeros((len(self.y1)-1))
        self.ymean=np.zeros((len(self.y1)-1))
        self.Lin=np.zeros((len(self.y1)-1))
        self.Lout=np.zeros((len(self.y1)-1))
        for i in range(len(self.x1)-1):
            self.vx[i,:]=(self.x1[i+1,:]-self.x1[i,:])/self.dt
            self.vy[i,:]=(self.y1[i+1,:]-self.y1[i,:])/self.dt
            self.xmean[i]=self.x1[i,:].mean()
            self.ymean[i]=self.y1[i,:].mean()
            self.Lin[i]=abs(((self.x1[i,:]-self.xmean[i])*self.vy[i,:]-(self.y1[i,:]-self.ymean[i])*self.vx[i,:]).mean())

    def getMSD(self):
        self.lenge=[]    
        for i in range(1,len(self.x1)-1):
            self.lenge.append(((self.x1[i,:]**2+self.y1[i,:]**2+self.x2[i,:]**2+self.y2[i,:]**2+self.x3[i,:]**2+self.y3[i,:]**2+self.x4[i,:]**2+self.y4[i,:]**2).mean()-(self.x1[i,:].mean()**2+self.y1[i,:].mean()**2+self.x2[i,:].mean()**2+self.y2[i,:].mean()**2+self.x3[i,:].mean()**2+self.y3[i,:].mean()**2+self.x4[i,:].mean()**2+self.y4[i,:].mean()**2))/4)#集群的转动惯量
        
        plt.plot(np.arange(1,len(self.lenge)+1)*self.dt,self.lenge,label="T:"+self.strn)
        plt.yscale('log')
        plt.xscale('log')
        #plt.legend(loc="upper left")
        plt.xlabel("$t/τ$")
        plt.ylabel("$MSD$")
        plt.title("$MSD$")
        plt.show()
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
        plt.show()
        #plt.legend(loc="upper left")
        self.dotmean=np.array(dot)[-200:].mean()
        return self.dotmean
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
        plt.show()
        #plt.legend(loc="upper left")
        self.cross=cross
        self.crossmean=np.array(cross)[-200:].mean()
        return self.crossmean
    def getFig(self):
        for i in range(0,20,1):
    #        for i in range(0,10000,1000):
            for j in range(len(self.x1[0,:])):
                plt.plot([self.x1[i,j],self.x2[i,j],self.x3[i,j],self.x4[i,j]],[self.y1[i,j],self.y2[i,j],self.y3[i,j],self.y4[i,j]],lw=8,solid_capstyle="round",antialiased=True)
                plt.scatter(self.x1[i,j],self.y1[i,j],s=100)
            plt.xlim(self.x1[i,:].min()-20,self.x1[i,:].max()+20)
            plt.ylim(self.y1[i,:].min()-15,self.y1[i,:].max()+15)
            plt.show()
    def getCrow(self):
        MSD=(np.log10(self.lenge))
        crow=np.where(MSD<2*MSD.min())[0]
        tau=(crow.max()-crow.min())*self.dt
        return MSD.min(),tau
    def getTrace(self):
        xmean=list(self.xmean)
        ymean=list(self.ymean)
        tar=self.getTraceR()
        for i in range(len(self.x1)-1):
            self.Lout[i]=abs(((self.x1[i,:]-self.rx)*self.vy[i,:]-(self.y1[i,:]-self.ry)*self.vx[i,:]).mean())
        plt.scatter(tar[0],tar[1],s=20,c="red")   
        plt.scatter(xmean,ymean,s=1)    
    def getTraceR(self):
        #print(x)
        tar=np.array([0.1,0.1,0.1])
        #tarUpdate=np.array([0,0,0])
        alpha=0.0005
        for i in range(2000):
            d=((self.xmean[10:]-tar[0])**2+(self.ymean[10:]-tar[1])**2-tar[2]**2)
            tar-=alpha*np.array([-2*((self.xmean[10:]-tar[0])*d).mean(),(-2*(self.ymean[10:]-tar[1])*d).mean(),-3*(2*d.mean())*tar[2]])
        self.rx=tar[0]
        self.ry=tar[1]
        self.r=tar[2]
        return tar
    def isBreak(self):
        if self.lenge[-1]>self.lenge[1]:
            return True
        else:
            return False
    def getBreakT(self):
        return len(np.where(self.lenge<2*min(self.lenge))[0])
        #return 0    
    def getDBSCAN(self,tar):
        dbscan=DBSCAN(4*3,min_samples=12)
        L=[]
        N=[]
        W=[]
        cross=[]
        dot=[]
        xmean=[]
        ymean=[]
        #tar=1
        for i in range(len(self.x1)*1//10,len(self.x1)-1,1):
            try:
                inpt=np.array([self.x1[i,:],self.y1[i,:],self.x2[i,:],self.y2[i,:],self.x3[i,:],self.y3[i,:],self.x4[i,:],self.y4[i,:]]).T
                dbscan.fit(inpt)
                labels=dbscan.labels_
                
                a=np.array(labels)   
                unique, counts = np.unique(a, return_counts=True)
                di=dict(zip(counts,unique))
                #tar=di[0]  if di[0]!=-1 else di[1]
                so=np.sort(counts)
                tar=di[so[-1]] if di[so[-1]]!=-1 else di[so[-2]]
                print(tar)
    
                TAR=np.where(labels==tar)[0]
                
                xmean.append(self.x1[i,TAR].mean())
                ymean.append(self.y1[i,TAR].mean())
                x=self.x1[i,TAR]-xmean[-1]
                y=self.y1[i,TAR]-ymean[-1]
                xo=self.x1[i-1,TAR]-xmean[-2]
                yo=self.y1[i-1,TAR]-ymean[-2]
                
                tempx=x#/np.sqrt(x**2+y**2)
                tempy=y#/np.sqrt(x**2+y**2)
                w=np.arccos((x*xo+y*yo)/np.sqrt(x**2+y**2)/np.sqrt(xo**2+yo**2)).mean()/self.dt

                P=np.array([np.array(self.x1[i,TAR]-self.x4[i,TAR]),np.array(self.y1[i,TAR]-self.y4[i,TAR])])
                R=np.array([np.array((self.x4[i,TAR]+self.x1[i,TAR])/2-self.x1[i,TAR].mean()),np.array((self.y4[i,TAR]+self.y1[i,TAR])/2-self.y1[i,TAR].mean())])
                Pmod=np.sqrt(P[0,:]**2+P[1,:]**2)
                Rmod=np.sqrt(R[0,:]**2+R[1,:]**2)
                P=P/Pmod
                R=R/Rmod
                per=[]
                for j in range(len(P[0,:])):
                    per.append(np.abs(np.linalg.det(np.array([P[:,j],R[:,j]]))))
                dot.append((((P*R)[0,:]+(P*R)[1,:])).mean())
                cross.append((np.array(per)).mean())
    
                L.append((((tempx)*self.vy[i,TAR]-(tempy)*self.vx[i,TAR]).mean()))
                N.append(len(TAR))
                W.append(w)
                #W.append(np.arccos((x*xo+y*yo)/np.sqrt(x**2+y**2)/np.sqrt(xo**2+yo**2)).mean()/self.dt)
                #L.append((((self.x1[i,np.where(labels==tar)]-xmean[-1])*self.vy[i,np.where(labels==tar)]-(self.y1[i,np.where(labels==tar)]-ymean[-1])*self.vx[i,np.where(labels==tar)]).mean()))
                
            except:
                pass
            
            
        #print("i:"+str(i))
        i=-1
        plt.scatter(self.x1[i,:],self.y1[i,:],s=20,c='red')
        plt.scatter(self.x2[i,:],self.y2[i,:],s=20,c='green')
        plt.scatter(self.x3[i,:],self.y3[i,:],s=20,c='green')
        plt.scatter(self.x4[i,:],self.y4[i,:],s=20,c='green')
               
        
        
        plt.scatter(self.x1[i,np.where(labels==tar)],self.y1[i,np.where(labels==tar)],s=20,c='red')
        plt.scatter(self.x2[i,np.where(labels==tar)],self.y2[i,np.where(labels==tar)],s=20,c='yellow')
        plt.scatter(self.x3[i,np.where(labels==tar)],self.y3[i,np.where(labels==tar)],s=20,c='yellow')
        plt.scatter(self.x4[i,np.where(labels==tar)],self.y4[i,np.where(labels==tar)],s=20,c='yellow')
        for j in range(len(self.x1[0,:])):
            plt.plot([self.x1[i,j],self.x2[i,j],self.x3[i,j],self.x4[i,j]],[self.y1[i,j],self.y2[i,j],self.y3[i,j],self.y4[i,j]] )
        plt.xlim(-10,110)
        plt.ylim(-10,110)
        plt.savefig("ana/pattern"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        plt.plot(L)
        #plt.ylim(-5,5)
        plt.savefig("ana/L"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        plt.plot(N)
        plt.savefig("ana/N"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        plt.plot(W,linewidth=1)
        plt.plot(np.zeros((len(W))))
        plt.savefig("ana/W"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        plt.scatter(N,W,s=0.1)
        plt.savefig("ana/WN"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        plt.scatter(N,L,s=0.1)
        plt.savefig("ana/LN"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        plt.scatter(dot,W,s=0.1)
        plt.savefig("ana/DW"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        plt.scatter(cross,W,s=0.1)
        plt.savefig("ana/CW"+self.stri+"/T"+self.strn+"Ts"+self.strm+"v0"+self.strv+".jpg")
        plt.close()
        W=np.array(W)
        t=np.where(~np.isnan(W))[0]
        return W[t][-100:].mean()
    def run(self):
        tem[self.id]=self.getDBSCAN(0)
        


w=[]
x=[]

mul=1
numb=48
for i in range(numb):
    tem=np.zeros(mul)#重复数
    task=[]
    for j in range(mul):
        task.append(read(j,0,i*0.02,10))
        task[-1].start()
    for t in task:
        t.join()
    w.append(np.mean(tem))
    x.append(i*0.02)
plt.plot(x,w)
plt.savefig("W.jpg")
np.savetxt("./W.txt",w)

