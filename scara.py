#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 16 13:42:18 2018

@author: mahou
"""
from numpy import zeros,sin,cos,arctan2,concatenate,sqrt,newaxis,vectorize

class Robot():
    def __init__(self,Params,centroids):
        self.L1=Params[0]
        self.L2=Params[1]
        self.Hm=Params[2]
        self.wp=Params[3]
        self.w1=Params[4]
        self.w2=Params[5]
        self.I0=Params[6]
        self.I1=Params[7]
        self.R_0=centroids[0]
        self.R_1=centroids[1]
        self.R_2=centroids[2]
        self.maxim=vectorize(self.maximo)
        pass
    def maximo(self,x):
        if(abs(x)<=1):
            return x
        else:
            return 1
    def fkine(self,q,centrid=False):
        n=q.shape[0]
        P0=zeros((n,3));P1=zeros((n,3));P2=zeros((n,3))
        P3=zeros((n,3))
        if(centrid):
            r1=sqrt(self.R_1[0]**2+self.R_1[1]**2)
            r2=sqrt(self.R_2[0]**2+self.R_2[1]**2)
            df1=arctan2(self.R_1[1],self.R_1[0])
            df2=arctan2(self.R_2[1],self.R_2[0])
            
            P0[:,0]=self.R_0[0]
            P0[:,1]=self.R_0[1]
            P0[:,2]=q[:,0]+self.R_0[2]
        
            P1[:,0]=r1*cos(q[:,1]+df1)
            P1[:,1]=r2*sin(q[:,1]+df1)
            P1[:,2]=q[:,0]+self.R_1[2]
            
            P3[:,0]=r2*cos(q[:,1]+q[:,2]+df2)
            P3[:,1]=r2*sin(q[:,1]+q[:,2]+df2)
            P3[:,2]=self.R_2[2]
            
            P2[:,0]=self.L1*cos(q[:,1])+P3[:,0]
            P2[:,1]=self.L1*sin(q[:,1])+P3[:,1]
            P2[:,2]=q[:,0]+P3[:,2]
        else:
            P0[:,2]=q[:,0]
        
            P1[:,0]=self.L1*cos(q[:,1])
            P1[:,1]=self.L1*sin(q[:,1])
            P1[:,2]=P0[:,2]
            
            P3[:,0]=self.L2*cos(q[:,1]+q[:,2])
            P3[:,1]=self.L2*sin(q[:,1]+q[:,2])
            P3[:,2]=0
            
            P2[:,0]=P1[:,0]+P3[:,0]
            P2[:,1]=P1[:,1]+P3[:,1]
            P2[:,2]=P0[:,2]+P3[:,2]
            pass
        return {'P0':P0,'P1':P1,'P2':P2,'P3':P3}
    def ikine(self,P):
        cs=(P[:,0]**2+P[:,1]**2-self.L1**2-self.L2**2)/(2*self.L1*self.L2)
        th=1-cs**2
        sn=sqrt(th)
        q2=arctan2(sn,cs)
        k1=self.L1+self.L2*cs
        k2=self.L2*sn
        q1=arctan2(P[:,1],P[:,0])-arctan2(k2,k1)
        return concatenate((P[:,2][:,newaxis],q1[:,newaxis],q2[:,newaxis]),axis=1)
    def jac(self,P):
        n=P['P1'].shape[0]
        J0=zeros((n,3,1))
        J1=zeros((n,3,2))
        J2=zeros((n,3,3))
        tmp=zeros((n,3))
        J0[:,2,:]=1
        
        J1[:,:,0]=J0[:,:,0]
        tmp[:,0]=-P['P1'][:,1]
        tmp[:,1]=P['P1'][:,0]
        J1[:,:,1]=tmp
        
        J2[:,:,0]=J0[:,:,0]
        tmp=0*tmp
        tmp[:,0]=-P['P2'][:,1]
        tmp[:,1]=P['P2'][:,0]
        J2[:,:,1]=tmp
        tmp=0*tmp
        tmp[:,0]=-P['P3'][:,1]
        tmp[:,1]=P['P3'][:,0]
        J2[:,:,2]=tmp
        return {'J0':J0,'J1':J1,'J2':J2}
    def djac(self,Qp,P,centroid=False):
        n=Qp.shape[0]
        J1=zeros((n,3,2))
        J2=zeros((n,3,3))
        tmp=zeros((n,3))
        
        if(centroid):
            tmp[:,0]=-P['P1c'][:,0]*Qp[:,1]
            tmp[:,1]=-P['P1c'][:,1]*Qp[:,1]
            J1[:,:,1]=tmp
            
            tmp=0*tmp
            tmp[:,0]=-P['P1'][:,0]*Qp[:,1]-P['P3'][:,0]*(Qp[:,1]+Qp[:,2])
            tmp[:,1]=-P['P1'][:,1]*Qp[:,1]-P['P3'][:,1]*(Qp[:,1]+Qp[:,2])
            J2[:,:,1]=tmp
            tmp=0*tmp
            tmp[:,0]=-P['P3'][:,0]*(Qp[:,1]+Qp[:,2])
            tmp[:,1]=-P['P3'][:,1]*(Qp[:,1]+Qp[:,2])
            J2[:,:,2]=tmp
        else:
            tmp[:,0]=-P['P1'][:,0]*Qp[:,1]
            tmp[:,1]=-P['P1'][:,1]*Qp[:,1]
            J1[:,:,1]=tmp
            
            tmp=0*tmp
            tmp[:,0]=-P['P1'][:,0]*Qp[:,1]-P['P3'][:,0]*(Qp[:,1]+Qp[:,2])
            tmp[:,1]=-P['P1'][:,1]*Qp[:,1]-P['P3'][:,1]*(Qp[:,1]+Qp[:,2])
            J2[:,:,1]=tmp
            tmp=0*tmp
            tmp[:,0]=-P['P3'][:,0]*(Qp[:,1]+Qp[:,2])
            tmp[:,1]=-P['P3'][:,1]*(Qp[:,1]+Qp[:,2])
            J2[:,:,2]=tmp
        return {'J1':J1,'J2':J2}
        
    pass


