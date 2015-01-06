
from __future__ import division
import numpy as np
from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from PIL import Image
import PngImagePlugin
import ttk
import sys
import os
from fractions import Fraction
from tkFileDialog import *
#import pdb

pi=np.pi


#Image._initialized=2
###################################################################"
##### Fonction projection  sur l'abaque
####################################################################

def proj(x,y,z): 
  
    if z==1: 
        X=0
        Y=0
    elif z<-0.000001:
        X=250
        Y=250
    else: 
            
        X=x/(1+z)
        Y=y/(1+z)
    
    return np.array([X,Y],float) 
    
###################################################################"
##### Fonction rotation 
####################################################################

def rotation(phi1,phi,phi2):
   phi1=phi1*pi/180;
   phi=phi*pi/180;
   phi2=phi2*pi/180;
   R=np.array([[np.cos(phi1)*np.cos(phi2)-np.cos(phi)*np.sin(phi1)*np.sin(phi2),
            -np.cos(phi)*np.cos(phi2)*np.sin(phi1)-np.cos(phi1)*
            np.sin(phi2),np.sin(phi)*np.sin(phi1)],[np.cos(phi2)*np.sin(phi1)
            +np.cos(phi)*np.cos(phi1)*np.sin(phi2),np.cos(phi)*np.cos(phi1)
            *np.cos(phi2)-np.sin(phi1)*np.sin(phi2), -np.cos(phi1)*np.sin(phi)],
            [np.sin(phi)*np.sin(phi2), np.cos(phi2)*np.sin(phi), np.cos(phi)]],float)
   return R

####################################################################
##### Fonction rotation autour d'un axe 
####################################################################

def Rot(th,a,b,c):
   th=th*pi/180;
   aa=a/np.linalg.norm([a,b,c]);
   bb=b/np.linalg.norm([a,b,c]);
   cc=c/np.linalg.norm([a,b,c]);
   c1=np.array([[1,0,0],[0,1,0],[0,0,1]],float)
   c2=np.array([[aa**2,aa*bb,aa*cc],[bb*aa,bb**2,bb*cc],[cc*aa,
                cc*bb,cc**2]],float)
   c3=np.array([[0,-cc,bb],[cc,0,-aa],[-bb,aa,0]],float)
   R=np.cos(th)*c1+(1-np.cos(th))*c2+np.sin(th)*c3

   return R    


####################################################################
##### Fonction cristal
####################################################################
def crist():
    global axesA,axeshA,axesB,axeshB,D,Dstar,V
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    e=eval(e_entry.get())
    d2=eval(d_label_var.get())
    alp=alp*pi/180;
    bet=bet*pi/180;
    gam=gam*pi/180;
    V=a*b*c*np.sqrt(1-(np.cos(alp)**2)-(np.cos(bet))**2-(np.cos(gam))**2+2*b*c*np.cos(alp)*np.cos(bet)*np.cos(gam))
    D=np.array([[a,b*np.cos(gam),c*np.cos(bet)],[0,b*np.sin(gam),  c*(np.cos(alp)-np.cos(bet)*np.cos(gam))/np.sin(gam)],[0,0,V/(a*b*np.sin(gam))]])
    Dstar=np.transpose(np.linalg.inv(D))
    G=np.array([[a**2,a*b*np.cos(gam),a*c*np.cos(bet)],[a*b*np.cos(gam),b**2,b*c*np.cos(alp)],[a*c*np.cos(bet),b*c*np.cos(alp),c**2]])    
    axes=np.zeros(((2*e+1)**3-1,3))
    axesh=np.zeros(((2*e+1)**3-1,3))
    id=0
    for i in range(-e,e+1):
        for j in range(-e,e+1):
            for k in range(-e,e+1):
                if (i,j,k)!=(0,0,0):
                    d=1/(np.sqrt(np.dot(np.array([i,j,k]),np.dot(np.linalg.inv(G),np.array([i,j,k])))))
                    if d>d2*0.1*np.amax([a,b,c]):
                        if var_uvw.get()==0:                    
                            axesh[id,:]=np.dot(Dstar,np.array([i,j,k],float))
                            axes[id,:]=np.array([i,j,k],float)
                        else:
                            axesh[id,:]=np.dot(D,np.array([i,j,k],float))
                            axes[id,:]=np.array([i,j,k],float)
                        id=id+1
    axesA=axes
    axesB=axes
    axeshA=axesh
    axeshB=axesh               
    return axesA,axeshA,axesB,axeshB,D,Dstar,V


def dm():
    global dmip
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)      
    dmip=dmip-eval(d_entry.get())
    d_label_var.set(dmip)
    crist()
    trace()
    
    
    return dmip
def dp():
    global dmip
    a=f.add_subplot(111)    
    a.figure.clear()
    a=f.add_subplot(111)      
    dmip=dmip+eval(d_entry.get())
    d_label_var.set(dmip)    
    crist()    
    trace()
   
    
    return dmip 

    

####################################################################
##### Fonction ajouter un pole
####################################################################
def poleA(pole1,pole2,pole3):
    global MA,axesA,axeshA,Ta,V,D,Dstar
    
    fp=f.add_subplot(111)
    Gs=np.array([pole1,pole2,pole3],float)
       
    Pp=np.zeros((1,2),float)
    
    if var_uvw.get()==0:                    
        Gsh=np.dot(Dstar,Gs)/np.linalg.norm(np.dot(Dstar,Gs))
    else:
        Gsh=np.dot(D,Gs)/np.linalg.norm(np.dot(D,Gs))
     
    S=np.dot(MA,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1=-pole1
        pole2=-pole2
        pole3=-pole3
    Pp=proj(S[0],S[1],S[2])*600/2
    l=str(int(pole1))+str(int(pole2))+str(int(pole3))
    fp.plot(Pp[0]+600/2,Pp[1]+600/2,'ro')    
    fp.annotate(l,(Pp[0]+600/2,Pp[1]+600/2))
    fp.axis([0,600,0,600])
    fp.axis('off')   
    fp.figure.canvas.draw() 
    
    axesA=np.vstack((axesA,np.array([pole1,pole2,pole3])))
    axesA=np.vstack((axesA,np.array([-pole1,-pole2,-pole3])))
    Ta=np.vstack((Ta,np.array([S[0],S[1],S[2]])))
    Ta=np.vstack((Ta,np.array([-S[0],-S[1],-S[2]])))
    axeshA=np.vstack((axeshA,np.array([Gsh[0],Gsh[1],Gsh[2]])))
    axeshA=np.vstack((axeshA,np.array([-Gsh[0],-Gsh[1],-Gsh[2]])))
    return axesA,axeshA,Ta
    
def poleB(pole1,pole2,pole3):
    global MB,axesB,axeshB,Tb,V,D,Dstar
    
    fp=f.add_subplot(111)
    Gs=np.array([pole1,pole2,pole3],float)
       
    Pp=np.zeros((1,2),float)
    
    if var_uvw.get()==0:                    
        Gsh=np.dot(Dstar,Gs)/np.linalg.norm(np.dot(Dstar,Gs))
    else:
        Gsh=np.dot(D,Gs)/np.linalg.norm(np.dot(D,Gs))
     
    S=np.dot(MB,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1=-pole1
        pole2=-pole2
        pole3=-pole3
    Pp=proj(S[0],S[1],S[2])*600/2
    l=str(int(pole1))+str(int(pole2))+str(int(pole3))
    fp.plot(Pp[0]+600/2,Pp[1]+600/2,'ro')    
    fp.annotate(l,(Pp[0]+600/2,Pp[1]+600/2))
    fp.axis([0,600,0,600])
    fp.axis('off')   
    fp.figure.canvas.draw() 
    
    axesB=np.vstack((axesB,np.array([pole1,pole2,pole3])))
    axesB=np.vstack((axesB,np.array([-pole1,-pole2,-pole3])))
    Tb=np.vstack((Tb,np.array([S[0],S[1],S[2]])))
    Tb=np.vstack((Tb,np.array([-S[0],-S[1],-S[2]])))
    axeshB=np.vstack((axeshB,np.array([Gsh[0],Gsh[1],Gsh[2]])))
    axeshB=np.vstack((axeshB,np.array([-Gsh[0],-Gsh[1],-Gsh[2]])))
    return axesB,axeshB,Tb
    
    
def addpoleA_sym():
        
    pole1A=eval(pole1A_entry.get())
    pole2A=eval(pole2A_entry.get())
    pole3A=eval(pole3A_entry.get())
    poleA(pole1A,pole2A,pole3A)
    poleA(pole1A,pole2A,-pole3A)
    poleA(pole1A,-pole2A,pole3A)
    poleA(-pole1A,pole2A,pole3A)
    poleA(pole2A,pole1A,pole3A)
    poleA(pole2A,pole1A,-pole3A)
    poleA(pole2A,-pole1A,pole3A)
    poleA(-pole2A,pole1A,pole3A)
    poleA(pole2A,pole3A,pole1A)
    poleA(pole2A,pole3A,-pole1A)
    poleA(pole2A,-pole3A,pole1A)
    poleA(-pole2A,pole3A,pole1A)
    poleA(pole1A,pole3A,pole2A)
    poleA(pole1A,pole3A,-pole2A)
    poleA(pole1A,-pole3A,pole2A)
    poleA(-pole1A,pole3A,pole2A)
    poleA(pole3A,pole1A,pole2A)
    poleA(pole3A,pole1A,-pole2A)
    poleA(pole3A,-pole1A,pole2A)
    poleA(-pole3A,pole1A,pole2A)
    poleA(pole3A,pole2A,pole1A)
    poleA(pole3A,pole2A,-pole1A)
    poleA(pole3A,-pole2A,pole1A)
    poleA(-pole3A,pole2A,pole1A)
    trace()

def addpoleB_sym():
        
    pole1B=eval(pole1B_entry.get())
    pole2B=eval(pole2B_entry.get())
    pole3B=eval(pole3B_entry.get())
    poleB(pole1B,pole2B,pole3B)
    poleB(pole1B,pole2B,-pole3B)
    poleB(pole1B,-pole2B,pole3B)
    poleB(-pole1B,pole2B,pole3B)
    poleB(pole2B,pole1B,pole3B)
    poleB(pole2B,pole1B,-pole3B)
    poleB(pole2B,-pole1B,pole3B)
    poleB(-pole2B,pole1B,pole3B)
    poleB(pole2B,pole3B,pole1B)
    poleB(pole2B,pole3B,-pole1B)
    poleB(pole2B,-pole3B,pole1B)
    poleB(-pole2B,pole3B,pole1B)
    poleB(pole1B,pole3B,pole2B)
    poleB(pole1B,pole3B,-pole2B)
    poleB(pole1B,-pole3B,pole2B)
    poleB(-pole1B,pole3B,pole2B)
    poleB(pole3B,pole1B,pole2B)
    poleB(pole3B,pole1B,-pole2B)
    poleB(pole3B,-pole1B,pole2B)
    poleB(-pole3B,pole1B,pole2B)
    poleB(pole3B,pole2B,pole1B)
    poleB(pole3B,pole2B,-pole1B)
    poleB(pole3B,-pole2B,pole1B)
    poleB(-pole3B,pole2B,pole1B)
    trace()
        
def addpoleA():
    
    pole1A=eval(pole1A_entry.get())
    pole2A=eval(pole2A_entry.get())
    pole3A=eval(pole3A_entry.get())
    poleA(pole1A,pole2A,pole3A)
    trace()
def addpoleB():
    
    pole1B=eval(pole1B_entry.get())
    pole2B=eval(pole2B_entry.get())
    pole3B=eval(pole3B_entry.get())
    poleB(pole1B,pole2B,pole3B)
    trace()
####################################################################
##### Fonction tracer plan
####################################################################
def trace_planA():
    global MA,axes,axesh,Ta,V,D,Dstar
    f2=f.add_subplot(111)    
    pole1A=eval(pole1A_entry.get())
    pole2A=eval(pole2A_entry.get())
    pole3A=eval(pole3A_entry.get())
    Gs=np.array([pole1A,pole2A,pole3A],float)
    if var_uvw.get()==0:                    
        Gsh=np.dot(Dstar,Gs)/np.linalg.norm(np.dot(Dstar,Gs))
    else:
        Gsh=np.dot(D,Gs)/np.linalg.norm(np.dot(D,Gs))
    S=np.dot(MA,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1A=-pole1A
        pole2A=-pole2A
        pole3A=-pole3A
    r=np.sqrt(S[0]**2+S[1]**2+S[2]**2)
    A=np.zeros((2,100))
    Q=np.zeros((1,2))
    if S[2]==0:
         t=90
         w=0
    else:
         t=np.arctan2(S[1],S[0])*180/pi
         w=0
    ph=np.arccos(S[2]/r)*180/pi
    for g in np.linspace(-pi,pi-0.00001,100):
        Aa=np.dot(Rot(t,0,0,1),np.dot(Rot(ph,0,1,0),np.array([np.sin(g),np.cos(g),0])))
        A[:,w]=proj(Aa[0],Aa[1],Aa[2])*600/2
        if A[0,w]<>75000:
            Q=np.vstack((Q,A[:,w]))
            w=w+1
    Q=np.delete(Q,0,0)    
    f2.plot(Q[:,0]+600/2,Q[:,1]+600/2,'r')
    f2.axis([0,600,0,600])
    f2.axis('off')      
    f2.figure.canvas.draw() 
    trace()

def trace_planB():
    global MB,axes,axesh,Tb,V,D,Dstar
    f2=f.add_subplot(111)    
    pole1B=eval(pole1B_entry.get())
    pole2B=eval(pole2B_entry.get())
    pole3B=eval(pole3B_entry.get())
    Gs=np.array([pole1B,pole2B,pole3B],float)
    if var_uvw.get()==0:                    
        Gsh=np.dot(Dstar,Gs)/np.linalg.norm(np.dot(Dstar,Gs))
    else:
        Gsh=np.dot(D,Gs)/np.linalg.norm(np.dot(D,Gs))
    S=np.dot(MB,Gsh)
    if S[2]<0:
        S=-S
        Gsh=-Gsh
        pole1B=-pole1B
        pole2B=-pole2B
        pole3B=-pole3B
    r=np.sqrt(S[0]**2+S[1]**2+S[2]**2)
    A=np.zeros((2,100))
    Q=np.zeros((1,2))
    if S[2]==0:
         t=90
         w=0
    else:
         t=np.arctan2(S[1],S[0])*180/pi
         w=0
    ph=np.arccos(S[2]/r)*180/pi
    for g in np.linspace(-pi,pi-0.00001,100):
        Aa=np.dot(Rot(t,0,0,1),np.dot(Rot(ph,0,1,0),np.array([np.sin(g),np.cos(g),0])))
        A[:,w]=proj(Aa[0],Aa[1],Aa[2])*600/2
        if A[0,w]<>75000:
            Q=np.vstack((Q,A[:,w]))
            w=w+1
    Q=np.delete(Q,0,0)    
    f2.plot(Q[:,0]+600/2,Q[:,1]+600/2,'r')
    f2.axis([0,600,0,600])
    f2.axis('off')      
    f2.figure.canvas.draw() 
    trace()

####################################################################
##### Click a pole
####################################################################    
def click_a_pole(event):
        
    global MB,Dstar,D
    x=event.x
    y=event.y
    x=(x-411)*2/620
    y=-(y-400)*2/620
    X=2*x/(1+x**2+y**2)
    Y=2*y/(1+x**2+y**2)
    Z=(-1+x**2+y**2)/(1+x**2+y**2)
    if Z<0:
        X=-X
        Y=-Y
    A=np.dot(np.linalg.inv(MB),np.array([X,Y,Z]))
    n=0
    L=np.zeros((3,16**3))                      
                           
                        
    for i in range(-8,9,1):
        for j in range(-8,9,1):
            for k in range(-8,9,1):
                if np.linalg.norm([i,j,k])<>0:
                    if var_uvw.get()==0:
                        Q=np.dot(Dstar,np.array([i,j,k],float))/np.linalg.norm(np.dot(Dstar,np.array([i,j,k],float)))
                        if np.abs(Q[0]-A[0])<0.05 and np.abs(Q[1]-A[1])<0.05 and np.abs(Q[2]-A[2])<0.05:
                            L[:,n]=np.array([i,j,k],float)
                            n=n+1
                           
                    else:
                          
                        Q=np.dot(D,np.array([i,j,k],float))/np.linalg.norm(np.dot(D,np.array([i,j,k],float)))
                        if np.abs(Q[0]-A[0])<0.05 and np.abs(Q[1]-A[1])<0.05 and np.abs(Q[2]-A[2])<0.05:
                            L[:,n]=np.array([i,j,k],float)
                            n=n+1
      

    if np.linalg.norm(L[:,0])<>0:
        poleB(L[0,0],L[1,0],L[2,0])
        trace()
####################################################################
##### Inclinaison-beta
####################################################################   
####################################################################
##### Fonction desorientation
####################################################################        
def Rota(t,u,v,w,g):
    Ae=np.dot(g,np.array([u,v,w]))
    Re=Rot(t,Ae[0],Ae[1],Ae[2])
    return Re

def cryststruct():
    global cs
    a=eval(a_entry.get())
    b=eval(b_entry.get())
    c=eval(c_entry.get())
    alp=eval(alp_entry.get())
    bet=eval(bet_entry.get())
    gam=eval(gam_entry.get())
    
    if  gam==90 and alp==90 and bet==90 and a==b and b==c:
       cs=1

    if gam==120 and alp==90 and bet==90:
        cs=2

    if gam==90 and alp==90 and bet==90 and a==b and b<>c: 
        cs=3  


    if alp<>90 and a==b and b==c:
        cs=4  

    if gam==90 and alp==90 and bet==90 and a<>b and b<>c:
        cs=5  
    
    if gam<>90 and alp==90 and bet==90 and a<>b and b<>c:
        cs=6  

    if gam<>90 and alp<>90 and bet<>90 and a<>b and b<>c: 
        cs=7  
    return cs
    
def Sy(g):
    global cs
    if cs==1:
        S1=Rota(90,1,0,0,g);
        S2=Rota(180,1,0,0,g);
        S3=Rota(270,1,0,0,g);
        S4=Rota(90,0,1,0,g);
        S5=Rota(180,0,1,0,g);
        S6=Rota(270,0,1,0,g);
        S7=Rota(90,0,0,1,g);
        S8=Rota(180,0,0,1,g);
        S9=Rota(270,0,0,1,g);
        S10=Rota(180,1,1,0,g);
        S11=Rota(180,1,0,1,g);
        
        S12=Rota(180,0,1,1,g);
        S13=Rota(180,-1,1,0,g);
        S14=Rota(180,-1,0,1,g);
        S15=Rota(180,0,-1,1,g);
        S16=Rota(120,1,1,1,g);
        S17=Rota(240,1,1,1,g);
        S18=Rota(120,-1,1,1,g);
        S19=Rota(240,-1,1,1,g);
        S20=Rota(120,1,-1,1,g);
        S21=Rota(240,1,-1,1,g);
        S22=Rota(120,1,1,-1,g);
        S23=Rota(240,1,1,-1,g);
        S24=np.eye(3,3);
        S=np.vstack((S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23,S24))
        
  
    
    if cs==2:
        S1=Rota(60,0,0,1,g);
        S2=Rota(120,0,0,1,g);
        S3=Rota(180,0,0,1,g);
        S4=Rota(240,0,0,1,g);
        S5=Rota(300,0,0,1,g);
        S6=np.eye(3,3);
        S7=Rota(180,0,0,1,g);
        S8=Rota(180,0,1,0,g);
        S9=Rota(180,1/2,np.sqrt(3)/2,0,g);
        S10=Rota(180,-1/2,np.sqrt(3)/2,0,g);
        S11=Rota(180,np.sqrt(3)/2,1/2,0,g);
        S12=Rota(180,-np.sqrt(3)/2,1/2,0,g);
        S=np.vstack((S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12))
        
       
        
    if cs==3:
        S1=Rota(90,0,0,1,g);
        S2=Rota(180,0,0,1,g);
        S3=Rota(270,0,0,1,g);
        S4=Rota(180,0,1,0,g);
        S5=Rota(180,1,0,0,g);
        S6=Rota(180,1,1,0,g);
        S7=Rota(180,1,-1,0,g);
        S8=np.eye(3,3)
        S=np.vstack((S1,S2,S3,S4,S5,S6,S7,S8))
        
        
      
        
    if cs==4:
        S1=Rota(60,0,0,1,g);
        S2=Rota(120,0,0,1,g);
        S3=Rota(180,0,0,1,g);
        S4=Rota(240,0,0,1,g);
        S5=Rota(300,0,0,1,g);
        S6=np.eye(3,3);
        S7=Rota(180,0,0,1,g);
        S8=Rota(180,0,1,0,g);
        S9=Rota(180,1/2,np.sqrt(3)/2,0,g);
        S10=Rota(180,-1/2,np.sqrt(3)/2,0,g);
        S11=Rota(180,np.sqrt(3)/2,1/2,0,g);
        S12=Rota(180,-np.sqrt(3)/2,1/2,0,g);
        S=np.vstack((S1,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12))
        
    
              
         
    if cs==5:
        S1=Rota(180,0,0,1,g);
        S2=Rota(180,1,0,0,g);
        S3=Rota(180,0,1,0,g);
        S4=np.eye(3,3);
        S=np.vstack((S1,S2,S3,S4))
        
                
        
    if cs==6:
        S1=Rota(180,0,1,0,g);
        S2=np.eye(3,3);
        S=np.vstack((S1,S2))
        
 
        
    if cs==7:
        S=np.eye(3,3);
        
    
    return S
    
def null(A, eps=1e-15):
    u, s, vh = np.linalg.svd(A)
    null_space = np.compress(s <= eps, vh, axis=0)
    return null_space.T
    
def desorientation():
    global D0,S,D1,cs,V,Qp
    a = f.add_subplot(111)
    a.figure.clear()
    a = f.add_subplot(111)
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    
    cryststruct()        
    phi1a=eval(phi1A_entry.get())
    phia=eval(phiA_entry.get())
    phi2a=eval(phi2A_entry.get())
    phi1b=eval(phi1B_entry.get())
    phib=eval(phiB_entry.get())
    phi2b=eval(phi2B_entry.get())
    
    gA=rotation(phi1a,phia,phi2a)
    gB=rotation(phi1b,phib,phi2b)
    k=0
    S=Sy(gA)
    
    D0=np.zeros((S.shape[0]/3,5))
    D1=np.zeros((S.shape[0]/3,3))
    Qp=np.zeros((S.shape[0]/3,2))
    for i in range(0,S.shape[0],3):
        In=np.dot(np.array([[S[i,0],S[i+1,0],S[i+2,0]],[S[i,1],S[i+1,1],S[i+2,1]],[S[i,2],S[i+1,2],S[i+2,2]]]),gA)
        Ing=np.dot(In,np.array([0,0,1]))
        In2=np.dot(Rot(-phi2b,Ing[0],Ing[1],Ing[2]),In)
        Ing2=np.dot(In2,np.array([1,0,0]))
        In3=np.dot(Rot(-phib,Ing2[0],Ing2[1],Ing2[2]),In2)
        Ing3=np.dot(In3,np.array([0,0,1]))
        A=np.dot(Rot(-phi1b,Ing3[0],Ing3[1],Ing3[2]),In3)-np.eye(3)
        V=null(A).T
        
        D0[k,0]=V[0,0]/np.linalg.norm(V)
        D0[k,1]=V[0,1]/np.linalg.norm(V)
        D0[k,2]=V[0,2]/np.linalg.norm(V)
       
        Ds1=np.dot(np.linalg.inv(gB),np.array([D0[k,0],D0[k,1],D0[k,2]]))
    
        F0=Fraction(Ds1[0]).limit_denominator(10)
        F1=Fraction(Ds1[1]).limit_denominator(10)
        F2=Fraction(Ds1[2]).limit_denominator(10)
                    
        D1[k,0]=F0.numerator*F1.denominator*F2.denominator
        D1[k,1]=F1.numerator*F0.denominator*F2.denominator
        D1[k,2]=F2.numerator*F0.denominator*F1.denominator
           
        D0[k,3]=np.arccos(0.5*(np.trace(A+np.eye(3))-1))*180/pi
        
        if D0[k,2]<0:
           D0[k,0]=-D0[k,0]
           D0[k,1]=-D0[k,1]
           D0[k,2]=-D0[k,2]
           D1[k,0]=-D1[k,0] 
           D1[k,1]=-D1[k,1]
           D1[k,2]=-D1[k,2]
           
       
       
        D0[k,4]=k
        Qp[k,:]=proj(D0[k,0],D0[k,1],D0[k,2])*600/2

        k=k+1
    
    a.plot(Qp[:,0]+600/2,Qp[:,1]+600/2,'ro')           
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw()
    trace()
    
    
    return Qp,S,D1
####################################################################
##### Fonction principale
####################################################################
def trace():
    global Ta,Tb,axesA,axeshA,MA,axesB,axeshB,MB,Qp,S,D1,show_ind,D0
    a = f.add_subplot(111)
    
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
        
    Pa=np.zeros((axesA.shape[0],2))
    Pb=np.zeros((axesB.shape[0],2))
    
    
    for i in range(0,axesA.shape[0]):
        axeshA[i,:]=axeshA[i,:]/np.linalg.norm(axeshA[i,:])
        Ta[i,:]=np.dot(MA,axeshA[i,:])
        Pa[i,:]=proj(Ta[i,0],Ta[i,1],Ta[i,2])*600/2
        if show_ind.get()==1:            
            m=np.amax([np.abs(axesA[i,0]),np.abs(axesA[i,1]),np.abs(axesA[i,2])])
            if (np.around(axesA[i,0]/m)==axesA[i,0]/m) & (np.around(axesA[i,1]/m)==axesA[i,1]/m) & (np.around(axesA[i,2]/m)==axesA[i,2]/m):
                sA=str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(int(axesA[i,2]/m))
            else:
               sA=str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(int(axesA[i,2]))  
            a.annotate(sA,(Pa[i,0]+600/2,Pa[i,1]+600/2))
    for i in range(0,axesB.shape[0]):
        axeshB[i,:]=axeshB[i,:]/np.linalg.norm(axeshB[i,:])
        Tb[i,:]=np.dot(MB,axeshB[i,:])
        Pb[i,:]=proj(Tb[i,0],Tb[i,1],Tb[i,2])*600/2
        if show_ind.get()==1:            
            m=np.amax([np.abs(axesB[i,0]),np.abs(axesB[i,1]),np.abs(axesB[i,2])])
            if (np.around(axesB[i,0]/m)==axesB[i,0]/m) & (np.around(axesB[i,1]/m)==axesB[i,1]/m) & (np.around(axesB[i,2]/m)==axesB[i,2]/m):
                sB=str(int(axesB[i,0]/m))+str(int(axesB[i,1]/m))+str(int(axesB[i,2]/m))
            else:
               sB=str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(int(axesB[i,2]))  
            a.annotate(sB,(Pb[i,0]+600/2,Pb[i,1]+600/2))   
    
    
    for l in range(0,int(S.shape[0]/3)):
        if show_angle.get()==1:
            sangle=str(np.round(D0[l,3],decimals=1))
            a.annotate(sangle,(Qp[l,0]+600/2,Qp[l,1]+600/2),size=8)
        if show_axe.get()==1:
            saxe=str(int(D1[l,0]))+','+str(int(D1[l,1]))+','+str(int(D1[l,2]))
            a.annotate(saxe,(Qp[l,0]+600/2,Qp[l,1]+600/2),size=8)
        if show_num.get()==1:
            snum=str(int(D0[l,4]))
            a.annotate(snum,(Qp[l,0]+600/2,Qp[l,1]+600/2),size=10)
        
        
        
    a.plot(Pa[:,0]+600/2,Pa[:,1]+600/2,'bo')           
    a.plot(Pb[:,0]+600/2,Pb[:,1]+600/2,'go')
    a.plot(Qp[:,0]+600/2,Qp[:,1]+600/2,'ro') 
    a.axis([0,600,0,600])
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw()
    

    
def princ():
    global Ta,Tb,MA,MB
    
    a = f.add_subplot(111)
    a.figure.clear()
    a = f.add_subplot(111)
    phi1a=eval(phi1A_entry.get())
    phia=eval(phiA_entry.get())
    phi2a=eval(phi2A_entry.get())
    phi1b=eval(phi1B_entry.get())
    phib=eval(phiB_entry.get())
    phi2b=eval(phi2B_entry.get())
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    crist()    
    Pa=np.zeros((axesA.shape[0],2))
    Ta=np.zeros((axesA.shape))
    Pb=np.zeros((axesB.shape[0],2))
    Tb=np.zeros((axesB.shape))
    
    for i in range(0,axesA.shape[0]):
        axeshA[i,:]=axeshA[i,:]/np.linalg.norm(axeshA[i,:])
        Ta[i,:]=np.dot(rotation(phi1a,phia,phi2a),axeshA[i,:])
        Pa[i,:]=proj(Ta[i,0],Ta[i,1],Ta[i,2])*600/2
        m=np.amax([np.abs(axesA[i,0]),np.abs(axesA[i,1]),np.abs(axesA[i,2])])
        if (np.around(axesA[i,0]/m)==axesA[i,0]/m) & (np.around(axesA[i,1]/m)==axesA[i,1]/m) & (np.around(axesA[i,2]/m)==axesA[i,2]/m):
            sA=str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(int(axesA[i,2]/m))
        else:
           sA=str(int(axesA[i,0]))+str(int(axesA[i,1]))+str(int(axesA[i,2]))  
        a.annotate(sA,(Pa[i,0]+600/2,Pa[i,1]+600/2))
    for i in range(0,axesB.shape[0]):
        axeshB[i,:]=axeshB[i,:]/np.linalg.norm(axeshB[i,:])
        Tb[i,:]=np.dot(rotation(phi1b,phib,phi2b),axeshB[i,:])
        Pb[i,:]=proj(Tb[i,0],Tb[i,1],Tb[i,2])*600/2
        m=np.amax([np.abs(axesB[i,0]),np.abs(axesB[i,1]),np.abs(axesB[i,2])])
        if (np.around(axesB[i,0]/m)==axesB[i,0]/m) & (np.around(axesB[i,1]/m)==axesB[i,1]/m) & (np.around(axesB[i,2]/m)==axesB[i,2]/m):
            sB=str(int(axesA[i,0]/m))+str(int(axesA[i,1]/m))+str(int(axesA[i,2]/m))
        else:
           sB=str(int(axesB[i,0]))+str(int(axesB[i,1]))+str(int(axesB[i,2]))  
        a.annotate(sB,(Pb[i,0]+600/2,Pb[i,1]+600/2))
                
               
    a.plot(Pa[:,0]+600/2,Pa[:,1]+600/2,'bo')
    a.plot(Pb[:,0]+600/2,Pb[:,1]+600/2,'go')
    a.axis([0,600,0,600])
    
    a.imshow(img)
    a.axis('off')   
    a.figure.canvas.draw() 
    
        
    MA=rotation(phi1a,phia,phi2a)
    MB=rotation(phi1b,phib,phi2b)
   
    return Ta,MA,MB,Tb
    
                    
######################################################################
# GUI
######################################################################

def file_save():
    global D1,D0
    fout = asksaveasfile(mode='w', defaultextension=".txt")
    
    for i in range(D1.shape[0]):
        text2save = str(int(D0[i,4]))+'\t'+'['+str(int(D1[i,0]))+','+str(int(D1[i,1]))+','+str(int(D1[i,2]))+']'+'\t '+str(np.around(D0[i,3],decimals=2))
        fout.write("%s\n" % text2save)
        
    fout.close()    

def image_save():
    
    s = asksaveasfile(mode='w', defaultextension=".jpg")    
    if s:    
        f.savefig(s.name)
    #s.close()
    
####################################################
#fonction d'initialisation
##################################################
def init():
    global var_uvw,D1,S,Qp,show_ind,show_angle,show_axe,show_num,dmip,d_label_var
    fn = os.path.join(os.path.dirname(__file__), 'stereo.png')      
    img=Image.open(fn)
    a = f.add_subplot(111)
    a.axis('off')
    a.imshow(img)
    a.figure.canvas.draw()
    S=np.zeros((1,5))
    Qp=np.zeros((1,2))
    D1=np.zeros((1,5))
    var_uvw=IntVar()
    show_ind=IntVar()
    show_angle=IntVar()
    show_axe=IntVar()
    show_num=IntVar()
    d_label_var=StringVar()
    d_label_var.set(0)
    dmip=0
    return var_uvw,show_ind,show_angle,show_axe,show_num
   

##############################################################
# fonction pour quitter
#######################################################
def _quit():
    root.quit()     # stops mainloop
    root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
#############################################################


root = Tk()
root.wm_title("Misorientation")
root.geometry('1220x798+10+40')
root.configure(bg = '#BDBDBD')
#root.resizable(0,0)
#s=ttk.Style()
#s.theme_use('clam')
style = ttk.Style()
theme = style.theme_use()
default = style.lookup(theme, 'background')



################################################
# Creation d'une zone pour tracer des graphiques
################################################
f = Figure(facecolor='white',figsize=[2,2],dpi=100)

canvas = FigureCanvasTkAgg(f, master=root)

canvas.get_tk_widget().place(x=0,y=0,height=800,width=800)
canvas._tkcanvas.bind('<Button-3>', click_a_pole)
canvas.show()
#toolbar = NavigationToolbar2TkAgg( canvas, root )
#toolbar.zoom('off')
#toolbar.update()


###################################################

init()
#import _imaging
#print _imaging.__file__

##############################################
# Boutons
##############################################
phi1A_entry = Entry (master=root)
phi1A_entry.place(relx=0.7,rely=0.5,relheight=0.03,relwidth=0.07)
phi1A_entry.configure(background="white")
phi1A_entry.configure(foreground="black")
phi1A_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phi1A_entry.configure(highlightcolor="#000000")
phi1A_entry.configure(insertbackground="#000000")
phi1A_entry.configure(selectbackground="#c4c4c4")
phi1A_entry.configure(selectforeground="black")

phiA_entry = Entry (master=root)
phiA_entry.place(relx=0.7,rely=0.55,relheight=0.03,relwidth=0.07)
phiA_entry.configure(background="white")
phiA_entry.configure(foreground="black")
phiA_entry.configure(highlightcolor="black")
phiA_entry.configure(insertbackground="black")
phiA_entry.configure(selectbackground="#c4c4c4")
phiA_entry.configure(selectforeground="black")

label_euler = Label (master=root)
label_euler.place(relx=0.77,rely=0.45,height=19,width=80)
label_euler.configure(activebackground="#cccccc")
label_euler.configure(activeforeground="black")
label_euler.configure(cursor="fleur")
label_euler.configure(foreground="black")
label_euler.configure(highlightcolor="black")
label_euler.configure(text='''Euler angles ''')

phi2A_entry = Entry (master=root)
phi2A_entry.place(relx=0.7,rely=0.6,relheight=0.03,relwidth=0.07)
phi2A_entry.configure(background="white")
phi2A_entry.configure(foreground="black")
phi2A_entry.configure(highlightcolor="black")
phi2A_entry.configure(insertbackground="black")
phi2A_entry.configure(selectbackground="#c4c4c4")
phi2A_entry.configure(selectforeground="black")

button_trace = Button (master=root)
button_trace.place(relx=0.7,rely=0.66,height=21,width=49)
button_trace.configure(activebackground="#f9f9f9")
button_trace.configure(activeforeground="black")
button_trace.configure(background="#ff0000")
button_trace.configure(command=princ)
button_trace.configure(foreground="black")
button_trace.configure(highlightcolor="black")
button_trace.configure(pady="0")
button_trace.configure(text='''Plot''')

Phi1A_label = Label (master=root)
Phi1A_label.place(relx=0.67,rely=0.5,height=19,width=33)
Phi1A_label.configure(activebackground="#cccccc")
Phi1A_label.configure(activeforeground="black")
Phi1A_label.configure(foreground="black")
Phi1A_label.configure(highlightcolor="black")
Phi1A_label.configure(text='''Phi1A''')

PhiA_label = Label (master=root)
PhiA_label.place(relx=0.67,rely=0.55,height=19,width=27)
PhiA_label.configure(activebackground="#cccccc")
PhiA_label.configure(activeforeground="black")
PhiA_label.configure(foreground="black")
PhiA_label.configure(highlightcolor="black")
PhiA_label.configure(text='''PhiA''')

Phi2A_label = Label (master=root)
Phi2A_label.place(relx=0.67,rely=0.6,height=19,width=33)
Phi2A_label.configure(activebackground="#cccccc")
Phi2A_label.configure(activeforeground="black")
Phi2A_label.configure(foreground="black")
Phi2A_label.configure(highlightcolor="black")
Phi2A_label.configure(text='''Phi2A''')

Cristal_label = Label (master=root)
Cristal_label.place(relx=0.66,rely=0.03,height=19,width=142)
Cristal_label.configure(text='''Crystal parameters''')

a_cristal_label = Label (master=root)
a_cristal_label.place(relx=0.68,rely=0.06,height=19,width=12)
a_cristal_label.configure(text='''a''')

b_cristal_label = Label (master=root)
b_cristal_label.place(relx=0.68,rely=0.1,height=19,width=12)
b_cristal_label.configure(activebackground="#f9f9f9")
b_cristal_label.configure(activeforeground="black")
b_cristal_label.configure(foreground="black")
b_cristal_label.configure(highlightcolor="black")
b_cristal_label.configure(text='''b''')

c_cristal_label = Label (master=root)
c_cristal_label.place(relx=0.68,rely=0.14,height=19,width=11)
c_cristal_label.configure(activebackground="#f9f9f9")
c_cristal_label.configure(activeforeground="black")
c_cristal_label.configure(foreground="black")
c_cristal_label.configure(highlightcolor="black")
c_cristal_label.configure(text='''c''')

alp_cristal_label = Label (master=root)
alp_cristal_label.place(relx=0.67,rely=0.19,height=19,width=32)
alp_cristal_label.configure(activebackground="#f9f9f9")
alp_cristal_label.configure(activeforeground="black")
alp_cristal_label.configure(foreground="black")
alp_cristal_label.configure(highlightcolor="black")
alp_cristal_label.configure(text='''alpha''')

bet_cristal_label = Label (master=root)
bet_cristal_label.place(relx=0.67,rely=0.23,height=19,width=28)
bet_cristal_label.configure(activebackground="#f9f9f9")
bet_cristal_label.configure(activeforeground="black")
bet_cristal_label.configure(foreground="black")
bet_cristal_label.configure(highlightcolor="black")
bet_cristal_label.configure(text='''beta''')

gam_cristal_label = Label (master=root)
gam_cristal_label.place(relx=0.66,rely=0.26,height=19,width=40)
gam_cristal_label.configure(activebackground="#f9f9f9")
gam_cristal_label.configure(activeforeground="black")
gam_cristal_label.configure(foreground="black")
gam_cristal_label.configure(highlightcolor="black")
gam_cristal_label.configure(text='''gamma''')

a_entry = Entry (master=root)
a_entry.place(relx=0.7,rely=0.06,relheight=0.03,relwidth=0.06)
a_entry.configure(background="white")
a_entry.configure(insertbackground="black")

b_entry = Entry (master=root)
b_entry.place(relx=0.7,rely=0.1,relheight=0.03,relwidth=0.06)
b_entry.configure(background="white")
b_entry.configure(foreground="black")
b_entry.configure(highlightcolor="black")
b_entry.configure(insertbackground="black")
b_entry.configure(selectbackground="#c4c4c4")
b_entry.configure(selectforeground="black")

c_entry = Entry (master=root)
c_entry.place(relx=0.7,rely=0.14,relheight=0.03,relwidth=0.06)
c_entry.configure(background="white")
c_entry.configure(foreground="black")
c_entry.configure(highlightcolor="black")
c_entry.configure(insertbackground="black")
c_entry.configure(selectbackground="#c4c4c4")
c_entry.configure(selectforeground="black")

alp_entry = Entry (master=root)
alp_entry.place(relx=0.7,rely=0.18,relheight=0.03,relwidth=0.06)
alp_entry.configure(background="white")
alp_entry.configure(foreground="black")
alp_entry.configure(highlightcolor="black")
alp_entry.configure(insertbackground="black")
alp_entry.configure(selectbackground="#c4c4c4")
alp_entry.configure(selectforeground="black")

bet_entry = Entry (master=root)
bet_entry.place(relx=0.7,rely=0.23,relheight=0.03,relwidth=0.06)
bet_entry.configure(background="white")
bet_entry.configure(foreground="black")
bet_entry.configure(highlightcolor="black")
bet_entry.configure(insertbackground="black")
bet_entry.configure(selectbackground="#c4c4c4")
bet_entry.configure(selectforeground="black")

gam_entry = Entry (master=root)
gam_entry.place(relx=0.7,rely=0.26,relheight=0.03,relwidth=0.06)
gam_entry.configure(background="white")
gam_entry.configure(foreground="black")
gam_entry.configure(highlightcolor="black")
gam_entry.configure(insertbackground="black")
gam_entry.configure(selectbackground="#c4c4c4")
gam_entry.configure(selectforeground="black")

uvw_button = Checkbutton (master=root)
uvw_button.place(relx=0.75,rely=0.66,relheight=0.03,relwidth=0.04)
uvw_button.configure(text='''uvw''')
uvw_button.configure(variable=var_uvw)

e_label = Label (master=root)
e_label.place(relx=0.66,rely=0.31,height=19,width=56)
e_label.configure(text='''max indice''')

e_entry = Entry (master=root)
e_entry.place(relx=0.72,rely=0.31,relheight=0.03,relwidth=0.05)
e_entry.configure(background="white")
e_entry.configure(insertbackground="black")


e2_label = Label (master=root)
e2_label.place(relx=0.68,rely=0.36,height=19,width=12)
e2_label.configure(text='''d''')

dm_button = Button (master=root)
dm_button.place(relx=0.7,rely=0.36,height=21,width=13)
dm_button.configure(activebackground="#f9f9f9")
dm_button.configure(activeforeground="black")
dm_button.configure(command=dm)
dm_button.configure(foreground="black")
dm_button.configure(highlightcolor="black")
dm_button.configure(pady="0")
dm_button.configure(text='''-''')

d_entry = Entry (master=root)
d_entry.place(relx=0.72,rely=0.36,relheight=0.02,relwidth=0.04)
d_entry.configure(background="white")
d_entry.configure(foreground="black")
d_entry.configure(highlightcolor="black")
d_entry.configure(insertbackground="black")
d_entry.configure(selectbackground="#c4c4c4")
d_entry.configure(selectforeground="black")

dp_button = Button (master=root)
dp_button.place(relx=0.76,rely=0.36,height=21,width=17)
dp_button.configure(activebackground="#f9f9f9")
dp_button.configure(activeforeground="black")
dp_button.configure(command=dp)
dp_button.configure(foreground="black")
dp_button.configure(highlightcolor="black")
dp_button.configure(pady="0")
dp_button.configure(text='''+''')

d_label = Label (master=root)
d_label.place(relx=0.73,rely=0.39,height=19,width=16)
d_label.configure(textvariable=d_label_var)


label_addpoleA = Label (master=root)
label_addpoleA.place(relx=0.81,rely=0.03,height=19,width=90)
label_addpoleA.configure(activebackground="#cccccc")
label_addpoleA.configure(activeforeground="black")
label_addpoleA.configure(foreground="black")
label_addpoleA.configure(highlightcolor="black")
label_addpoleA.configure(text='''Add pole A''')

pole1A_entry = Entry (master=root)
pole1A_entry.place(relx=0.81,rely=0.06,relheight=0.02
,relwidth=0.04)
pole1A_entry.configure(background="white")
pole1A_entry.configure(foreground="black")
pole1A_entry.configure(highlightcolor="black")
pole1A_entry.configure(insertbackground="black")
pole1A_entry.configure(selectbackground="#c4c4c4")
pole1A_entry.configure(selectforeground="black")

pole2A_entry = Entry (master=root)
pole2A_entry.place(relx=0.87,rely=0.06,relheight=0.02
,relwidth=0.04)
pole2A_entry.configure(background="white")
pole2A_entry.configure(foreground="black")
pole2A_entry.configure(highlightcolor="black")
pole2A_entry.configure(insertbackground="black")
pole2A_entry.configure(selectbackground="#c4c4c4")
pole2A_entry.configure(selectforeground="black")

pole3A_entry = Entry (master=root)
pole3A_entry.place(relx=0.93,rely=0.06,relheight=0.02
,relwidth=0.04)
pole3A_entry.configure(background="white")
pole3A_entry.configure(foreground="black")
pole3A_entry.configure(highlightcolor="black")
pole3A_entry.configure(insertbackground="black")
pole3A_entry.configure(selectbackground="#c4c4c4")
pole3A_entry.configure(selectforeground="black")

addpoleA_button = Button (master=root)
addpoleA_button.place(relx=0.81,rely=0.11,height=31,width=57)
addpoleA_button.configure(activebackground="#f9f9f9")
addpoleA_button.configure(activeforeground="black")
addpoleA_button.configure(command=addpoleA)
addpoleA_button.configure(foreground="black")
addpoleA_button.configure(highlightcolor="black")
addpoleA_button.configure(pady="0")
addpoleA_button.configure(text='''Add''')

symA_button = Button (master=root)
symA_button.place(relx=0.87,rely=0.11,height=31,width=71)
symA_button.configure(command=addpoleA_sym)
symA_button.configure(pady="0")
symA_button.configure(text='''Symmetry''')

trace_planA_button = Button (master=root)
trace_planA_button.place(relx=0.93,rely=0.11,height=31,width=61)
trace_planA_button.configure(command=trace_planA)
trace_planA_button.configure(pady="0")
trace_planA_button.configure(text='''Plane''')

phi1B_entry = Entry (master=root)
phi1B_entry.place(relx=0.84,rely=0.5,relheight=0.03,relwidth=0.07)
phi1B_entry.configure(background="white")

phi1B_entry.configure(foreground="black")
phi1B_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phi1B_entry.configure(highlightcolor="#000000")
phi1B_entry.configure(insertbackground="#000000")
phi1B_entry.configure(selectbackground="#c4c4c4")
phi1B_entry.configure(selectforeground="black")

Phi1B = Label (master=root)
Phi1B.place(relx=0.8,rely=0.5,height=29,width=33)
Phi1B.configure(activebackground="#cccccc")
Phi1B.configure(activeforeground="black")
Phi1B.configure(foreground="black")
Phi1B.configure(highlightcolor="black")
Phi1B.configure(text='''Phi1B''')

PhiB_label1 = Label (master=root)
PhiB_label1.place(relx=0.81,rely=0.55,height=19,width=26)
PhiB_label1.configure(activebackground="#cccccc")
PhiB_label1.configure(activeforeground="black")
PhiB_label1.configure(foreground="black")
PhiB_label1.configure(highlightcolor="black")
PhiB_label1.configure(text='''PhiB''')

Phi2B_label2 = Label (master=root)
Phi2B_label2.place(relx=0.8,rely=0.6,height=19,width=32)
Phi2B_label2.configure(activebackground="#cccccc")
Phi2B_label2.configure(activeforeground="black")
Phi2B_label2.configure(foreground="black")
Phi2B_label2.configure(highlightcolor="black")
Phi2B_label2.configure(text='''Phi2B''')

phiB_entry = Entry (master=root)
phiB_entry.place(relx=0.84,rely=0.55,relheight=0.03,relwidth=0.07)
phiB_entry.configure(background="white")

phiB_entry.configure(foreground="black")
phiB_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phiB_entry.configure(highlightcolor="#000000")
phiB_entry.configure(insertbackground="#000000")
phiB_entry.configure(selectbackground="#c4c4c4")
phiB_entry.configure(selectforeground="black")

phi2B_entry = Entry (master=root)
phi2B_entry.place(relx=0.84,rely=0.6,relheight=0.03,relwidth=0.07)
phi2B_entry.configure(background="white")

phi2B_entry.configure(foreground="black")
phi2B_entry.configure(highlightbackground="#e0e0dfdfe3e3")
phi2B_entry.configure(highlightcolor="#000000")
phi2B_entry.configure(insertbackground="#000000")
phi2B_entry.configure(selectbackground="#c4c4c4")
phi2B_entry.configure(selectforeground="black")

label_addpoleB = Label (master=root)
label_addpoleB.place(relx=0.81,rely=0.2,height=19,width=90)
label_addpoleB.configure(activebackground="#cccccc")
label_addpoleB.configure(activeforeground="black")
label_addpoleB.configure(foreground="black")
label_addpoleB.configure(highlightcolor="black")
label_addpoleB.configure(text='''Add pole B''')

pole1B_entry = Entry (master=root)
pole1B_entry.place(relx=0.81,rely=0.24,relheight=0.02
,relwidth=0.04)
pole1B_entry.configure(background="white")
pole1B_entry.configure(foreground="black")
pole1B_entry.configure(highlightcolor="black")
pole1B_entry.configure(insertbackground="black")
pole1B_entry.configure(selectbackground="#c4c4c4")
pole1B_entry.configure(selectforeground="black")

pole2B_entry = Entry (master=root)
pole2B_entry.place(relx=0.87,rely=0.24,relheight=0.02
,relwidth=0.04)
pole2B_entry.configure(background="white")
pole2B_entry.configure(foreground="black")
pole2B_entry.configure(highlightcolor="black")
pole2B_entry.configure(insertbackground="black")
pole2B_entry.configure(selectbackground="#c4c4c4")
pole2B_entry.configure(selectforeground="black")

pole3B_entry = Entry (master=root)
pole3B_entry.place(relx=0.93,rely=0.24,relheight=0.02
,relwidth=0.04)
pole3B_entry.configure(background="white")
pole3B_entry.configure(foreground="black")
pole3B_entry.configure(highlightcolor="black")
pole3B_entry.configure(insertbackground="black")
pole3B_entry.configure(selectbackground="#c4c4c4")
pole3B_entry.configure(selectforeground="black")

addpoleB_button = Button (master=root)
addpoleB_button.place(relx=0.81,rely=0.28,height=31,width=55)
addpoleB_button.configure(activebackground="#f9f9f9")
addpoleB_button.configure(activeforeground="black")
addpoleB_button.configure(command=addpoleB)
addpoleB_button.configure(foreground="black")
addpoleB_button.configure(highlightcolor="black")
addpoleB_button.configure(pady="0")
addpoleB_button.configure(text='''Add''')

symB_button = Button (master=root)
symB_button.place(relx=0.87,rely=0.28,height=31,width=70)
symB_button.configure(command=addpoleB_sym)
symB_button.configure(pady="0")
symB_button.configure(text='''Symmetry''')

trace_planB_button = Button (master=root)
trace_planB_button.place(relx=0.93,rely=0.28,height=31,width=59)
trace_planB_button.configure(command=trace_planB)
trace_planB_button.configure(pady="0")
trace_planB_button.configure(text='''Plane''')

button_desorientation = Button (master=root)
button_desorientation.place(relx=0.81,rely=0.66,height=21,width=98)

button_desorientation.configure(activebackground="#f9f9f9")
button_desorientation.configure(activeforeground="black")
button_desorientation.configure(background="#00ff00")
button_desorientation.configure(command=desorientation)
button_desorientation.configure(foreground="black")
button_desorientation.configure(highlightcolor="black")
button_desorientation.configure(pady="0")
button_desorientation.configure(text='''Misorientation''')

show_ind_button = Checkbutton (master=root)
show_ind_button.place(relx=0.81,rely=0.7,relheight=0.03
                ,relwidth=0.11)
show_ind_button.configure(text='''Show indices''')
show_ind_button.configure(variable=show_ind)

show_angle_button = Checkbutton (master=root)
show_angle_button.place(relx=0.81,rely=0.74,relheight=0.03
                ,relwidth=0.11)
show_angle_button.configure(text='''Show angle''')
show_angle_button.configure(variable=show_angle)

show_axe_button = Checkbutton (master=root)
show_axe_button.place(relx=0.81,rely=0.78,relheight=0.03
                ,relwidth=0.11)
show_axe_button.configure(text='''Show axes''')
show_axe_button.configure(variable=show_axe)

show_num_button = Checkbutton (master=root)
show_num_button.place(relx=0.81,rely=0.82,relheight=0.03
                ,relwidth=0.11)
show_num_button.configure(text='''Show Number''')
show_num_button.configure(variable=show_num)

menu = Menu(master=root)
filemenu = Menu(menu, tearoff=0)
menu.add_cascade(label="Save", menu=filemenu)

root.config(menu=menu)
filemenu.add_command(label="Save data", command=file_save) 
filemenu.add_command(label="Save figure", command=image_save) 
######################################################################################################
######## importer des structures cristallines depuis un fichier Nom,a,b,c,alpha,beta,gamma,space group
######################################################################################################

def structure(i0):
    global x0
    
    a_entry.delete(0,END)
    a_entry.insert(1,eval(x0[i0][1]))
    b_entry.delete(0,END)    
    b_entry.insert(1,eval(x0[i0][2]))
    c_entry.delete(0,END)    
    c_entry.insert(1,eval(x0[i0][3]))
    alp_entry.delete(0,END)    
    alp_entry.insert(1,eval(x0[i0][4]))
    bet_entry.delete(0,END)    
    bet_entry.insert(1,eval(x0[i0][5]))
    gam_entry.delete(0,END)    
    gam_entry.insert(1,eval(x0[i0][6]))
    
def createstructure(i):
    return lambda:structure(i)    
    
cristalmenu=Menu(menu,tearoff=0)
menu.add_cascade(label="Structures", menu=cristalmenu)
file_struct=open(os.path.join(os.path.dirname(__file__), 'structure.txt') ,"r")

x0=[]
i=0
for line in file_struct:
    x0.append(map(str, line.split()))
    cristalmenu.add_command(label=x0[i][0], command=createstructure(i))
    i=i+1
  
file_struct.close()

#######################################################################################################

phi1A_entry.insert(0,0)
phiA_entry.insert(0,0)
phi2A_entry.insert(0,0)
phi1B_entry.insert(0,0)
phiB_entry.insert(0,0)
phi2B_entry.insert(0,0)
e_entry.insert(1,1)
d_entry.insert(1,1)


mainloop()     
              
              
              
