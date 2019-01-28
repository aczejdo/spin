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


def markov_ising(S,E):
    J=S.shape[0]

    kx=np.random.randint(1,J-1)
    ky=np.random.randint(1,J-1)
    sigk = S[kx,ky]
    
    h=(S[kx+1,ky]+S[kx-1,ky]+S[kx,ky+1]+S[kx,ky-1])
    dE = 2*h*sigk
    Gamma= np.exp(-B*dE)
    """
    print('coord',kx,ky)
    print('E/n',E/N)
    print('G',Gamma)
    print('dE',dE)
    """
    if (np.random.rand() < Gamma):
        S[kx,ky] = -sigk
        E=E+dE
    
    return(S,E)

def init():
    #Tc=2/np.log(1+np.sqrt(2))
    T=2
    B=1/T	
    N=128
    maxit = int(1E5)
    
    Evec=np.zeros((maxit))
    
    M=np.zeros((N+2,N+2))
    M[1:N+1,1:N+1]=np.random.randint(0,2,size=(N,N))*2 - 1
    
    Mag=np.zeros((maxit))
    return(M,Evec,maxit,N,B,Mag)

def run_model(M,Energy,maxit,Evec):
    for i in range(maxit):
        M=bnc_periodic(M)
        M,Energy=markov_ising(M,Energy)
        Evec[i]=Energy
        Mag[i]=np.sum(M)
    return(M,Evec,Mag)


def plot_singlerun(Mag,Evec,M):
    fig = plt.figure()
    plt.imshow(M,cmap='bone')
    plt.show()    
    
    fig = plt.figure()
    plt.plot(np.arange(0,maxit),Evec)
    plt.title('energy')
    plt.show
    
    fig = plt.figure()
    plt.plot(np.arange(0,maxit),Mag)
    plt.title('magnetization')
    plt.show
    print('evec -1',Evec[-1],'\n')
    print('correct E',getE(M),'\n')
    print('evec/n',Evec[-1]/N,'getE/n',getE(M)/N,'\n')
    print('meanEvec/N',np.mean(Evec)/N,'\n')
    return()




M,Evec,maxit,N,B, Mag = init()
Energy=getE(M)
#======================
#plot initial config
#=======================
fig = plt.figure()
ax = fig.add_subplot(111)
plt.imshow(M,cmap='bone')
plt.show()



#=========================
#Run MC MH with bc
#==========================
M,Evec,Mag = run_model(M,Energy,maxit,Evec)


#=============================
#plot results
#==============================
plot_singlerun(Mag,Evec,M)
