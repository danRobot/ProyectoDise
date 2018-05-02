#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 19:52:53 2018

@author: mahou
"""

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy.interpolate import interp1d
from scipy import integrate

class Trayectoria():
    def __init__(self,points=0,time=0,Path=np.array([0,0,0])):
        self.n=points
        self.t=time
        self.P=Path
        _=self.measure()
        pass
    def measure(self,n=0):
        if(n==0):
            n=self.n
        self.L=-1
        x0=self.P[0,0]
        y0=self.P[0,1]
        z0=self.P[0,2]
        L=0
        for p in (self.P[1:n]):
            L=L+np.sqrt((p[0]-x0)**2+(p[1]-y0)**2+(p[2]-z0)**2)
            x0=p[0]
            y0=p[1]
            z0=p[2]
            pass
        self.L=L
        return self.L
    def measure2(self):
        self.L=-1
        x0=self.P[0,0]
        y0=self.P[0,1]
        z0=self.P[0,2]
        L=0
        i=0
        L1=np.zeros((self.n,))
        for p in (self.P):
            L=L+np.sqrt((p[0]-x0)**2+(p[1]-y0)**2+(p[2]-z0)**2)
            x0=p[0]
            y0=p[1]
            z0=p[2]
            L1[i]=L
            i=i+1
            pass
        self.L1=L1
        return L1
    def createVel(self,d):
        td=d*self.t
        self.time=np.linspace(0,self.t,self.n)
        #t1=self.time[0:int(self.n/2)]
        #t2=self.time[int(self.n/2):self.n]
        a=2/(self.t-td)
        tr=(self.t-td)
        ts=self.t-4*tr
        tr2=3*tr+ts
        """sec=int(self.n/7)
        self.time=np.concatenate([
                np.linspace(0,2*tr,3*sec),
                np.linspace(2*tr,0.75*tr2,sec),
                np.linspace(0.75*tr2,self.t,3*sec)])"""
        v1=1/(1+np.exp(-4*a*(self.time-tr)))
        v2=1-1/(1+np.exp(-4*a*(self.time-tr2)))
        v=v1*v2
        Vmax=self.L/np.trapz(v,dx=self.time[2]-self.time[1])
        self.v=Vmax*v
        return self.time,self.v
    def createVel2(self):
        Vmax=np.pi*self.t/self.L
        self.time=np.linspace(0,self.t,self.n)
        a=1-((self.time-0.5*self.t)**2)/self.t
        self.v=Vmax*np.sqrt(a)
        return self.time,self.v
    def assign(self):
        V=np.zeros((self.n,3))
        uv=np.zeros((self.n,3))
        l2=[]
        t2=[]
        self.x=integrate.cumtrapz(self.v,self.time,initial=0)
        fx=interp1d(self.x,self.time)
        fv=interp1d(self.time,self.v)
        for i in range(self.n):
            if(i==(self.n-1)):
                j=i
                norm=np.linalg.norm(self.P[j]-self.P[i-1])
                if(norm==0):
                    u=np.array([0,0,0])
                else:
                    u=(self.P[j]-self.P[i-1])/norm
            else:
                j=i+1
                norm=np.linalg.norm(self.P[j]-self.P[i])
                if(norm==0):
                    u=np.array([0,0,0])
                else:
                    u=(self.P[j]-self.P[i])/norm
                pass
            L0=self.L1[i]
            if((L0-self.x[self.n-1])>0):
                L0=self.x[self.n-1]
                pass
            ts=fx(L0)
            Vi=fv(ts)
            V[i,:]=u*Vi
            uv[i,:]=u
            l2.append(L0)
            t2.append(ts)
            pass
        return V,uv,np.array(t2),np.array(l2)
    pass
def createPath(Params,n,xi,xf,yi,yf):
    t=np.linspace(-2*pi,-0.5,n,dtype=np.float64)
    #X=np.abs(20+np.cos(t)+np.cos(2*t)+np.cos(4*t)+5*np.cos(6*t)+10*np.cos(8*t))
    #Y=(np.sin(t)+np.sin(3*t)+np.sin(5*t)+np.sin(7*t)+np.sin(9*t))
    #Z=np.abs(20+np.sin(t)+0.5*np.cos(2*t)+0.2*np.sin(3*t)+0.15*np.cos(5*t)+2*np.sin(7*t))
    X=np.linspace(xi,xf,n,dtype=np.float64)
    Y=np.linspace(yi,yf,n,dtype=np.float64)
    Z=np.linspace(-1,1,n,dtype=np.float64)
    #Z=-(Z-0.15)**2+0.15**2
    Z=0.25*np.sqrt(1**2-Z**4)
    #Z=1/(1+np.exp(-Z))
    #X=2*np.cos(t)+0.5
    #X=np.array(list(reversed(X)))
    #Y=2*np.sin(t)+0.5
    #Y=np.array(list(reversed(Y)))
    #Z=X*Y
    #ma=max(np.sqrt(X**2+Y**2+Z**2))*(np.sqrt(Params[0]**2+Params[1]**2+Params[2]**2))
    maz=max(Z)/Params[2]
    may=max(Y)/Params[1]
    ma=max(X)/Params[0]
    #X=X/ma
    #Y=Y/may
    #Z=Z/maz
    shape=(n,1)
    return np.concatenate([X.reshape(shape),Y.reshape(shape),Z.reshape(shape)],axis=1)

def createVel(d,tf,x,n):
    td=d*tf
    t=np.linspace(0,tf,n)
    a=2/(tf-td)
    tr=(tf-td)
    ts=tf-4*tr
    tr2=3*tr+ts
    v1=1/(1+np.exp(-4*a*(t-tr)))
    v2=1-1/(1+np.exp(-4*a*(t-tr2)))
    v=v1*v2
    Vmax=x/np.trapz(v,dx=t[2]-t[1])
    v=Vmax*v
    return t,v



def plot_colourline(fig,x,y,z,c):
    #c = cm.jet((c-np.min(c))/(np.max(c)-np.min(c)))
    ax = fig.gca(projection='3d')
    for i in np.arange(len(x)-1):
        ax.plot([x[i],x[i+1]], [y[i],y[i+1]],[z[i],z[i+1]], c=c[i])
    return

def measure(X,Y,Z):
    x0=X[0]
    y0=Y[0]
    z0=Z[0]
    L=0
    for x,y,z in zip(X[1:],Y[1:],Z[1:]):
        L=L+np.sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
        x0=x
        y0=y
        z0=z
        pass
    return L

def assign(P,Vmag):
    n=P.shape[0]
    V=np.zeros((n,3))
    uv=np.zeros((n,3))
    for i in range(n):
        if(i==(n-1)):
            j=i
            norm=np.linalg.norm(P[j]-P[i-1])
            if(norm==0):
                u=np.array([0,0,0])
            else:
                u=(P[j]-P[i-1])/norm
        else:
            j=i+1
            norm=np.linalg.norm(P[j]-P[i])
            if(norm==0):
                u=np.array([0,0,0])
            else:
                u=(P[j]-P[i])/norm
            pass
        V[i,:]=u*Vmag[i]
        uv[i,:]=u
        pass
    return V,uv

