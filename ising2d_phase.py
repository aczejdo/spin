# -*- coding: utf-8 -*-
"""
Created on Tue May 29 19:01:46 2018

@author: Aleks
"""
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

def bnc_periodic(A):
    #N+2 by N+2 in
    J=A.shape[0]-1
    A[0, 0:J] =  A[J-1, 0:J]
    A[J, 0:J] =  A[1  , 0:J]
    A[0:J, 0] =  A[0:J, J-1]
    A[0:J, J] =  A[0:J, 1]    
    return(A)
    

def getE(S):
    i = np.arange(1,S.shape[0]-1)
    j = np.arange(1,S.shape[1]-1)
    E=np.sum((S[i+1,j]+S[i-1,j]+S[i,j+1]+S[i,j-1])*S[i,j])
    J=1
    return(-J*E)


def markov_ising(S,E,B):
    J=S.shape[0]

    kx=np.random.randint(1,J-1)
    ky=np.random.randint(1,J-1)
    sigk = S[kx,ky]
    
    h=(S[kx+1,ky]+S[kx-1,ky]+S[kx,ky+1]+S[kx,ky-1])
    dE = 2*h*sigk
    Gamma= np.exp(-B*dE)

    if (np.random.rand() < Gamma):
        S[kx,ky] = -sigk
        E=E+dE
    
    return(S,E)

def init():
    #Tc=2/np.log(1+np.sqrt(2))
    N=32
    maxit = int(1E6)
    
    Evec=np.zeros((maxit))
    
    M=np.zeros((N+2,N+2))
    M[1:N+1,1:N+1]=np.random.randint(0,2,size=(N,N))*2 - 1
    
    Mag=np.zeros((maxit))
    return(M,Evec,maxit,N,Mag)

def run_model(M,Energy,maxit,Evec,B):
    for i in range(maxit):
        M=bnc_periodic(M)
        M,Energy=markov_ising(M,Energy,B)
        Evec[i]=Energy
        Mag[i]=np.sum(M)
    return(M,Evec,Mag)






#=========================
#Run MC MH with bc
#==========================
T=1E-7
runs=100
Tmax=6
dT=Tmax/runs
index=0


Magvec    =np.zeros((runs))
Tempvec   =np.arange(T,Tmax,dT)
Energyvec =np.zeros((runs))
while T<Tmax:
   
    B=1/Tempvec[index]
    M,Evec,maxit,N,Mag = init()
    Energy=getE(M)
    
    M,null1,Mag=run_model(M,Energy,maxit,Evec,B)
    Magvec[index] = np.abs(np.sum(M)/(N*N))
    Energyvec[index]=getE(M)
    
    T=T+dT
    print(index)
    index=index+1



fig =plt.figure()
plt.plot(Tempvec,Magvec)
plt.title('mamps')
plt.show()

fig =plt.figure()
plt.plot(Tempvec,Energyvec)
plt.title('mamps')
plt.show()
