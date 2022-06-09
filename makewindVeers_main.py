# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 18:35:06 2019

@author: radha
"""
from numpy import *
from scipy import *
from matplotlib import *

def ieckaimal(f, hubht, uhub, TI):
   Lambda = 0.7*min(60, hubht)    # TurbSim (23)
   L = 8.1*Lambda                 # TurbSim (22) 
   I = TI*(0.75 + 5.6/uhub)       # From IEC Ed 3 (Burton 2.23)
   sigma_u = I*uhub
   S = 4*pow(sigma_u,2)*L/uhub/pow((1 + 6*f*L/uhub),(5/3)) # TurbSim (21)
   return S

def makeKaimalSpectrum (hubht, uhub, TI_ref,dt_fft,Nf):   
    fk = logspace(log10(1/600),log10(0.5*1/dt_fft),Nf)    # selected frequencies
    G=ieckaimal(fk, hubht, uhub, TI_ref)      # PSD
    # get frequency bin width
    d_fft = diff(fk)/2
    edges = hstack((1e-10, fk[0:-1]+d_fft, fk[-1]+d_fft[-1]))
    dfk=diff(edges)
    var0=(TI_ref*(0.75 + 5.6/uhub)*uhub)**2     # specified variance
    s0=var0-sum(G*dfk)                        # missing variance
    Sres=ones(size(dfk))*s0/size(dfk)        # spread missing variance over all frequencies
    # the discrete spectrum (see notebook2015, p. 23)
    S0=G*dfk+Sres
    return S0,fk,edges

def makeKaimalSpectrum_630s (hubht, uhub, TI_ref,dt_fft,Nf):   
    fk = logspace(log10(1/630),log10(0.5*1/dt_fft),Nf)    # selected frequencies
    G=ieckaimal(fk, hubht, uhub, TI_ref)      # PSD
    # get frequency bin width
    d_fft = diff(fk)/2
    edges = hstack((1e-10, fk[0:-1]+d_fft, fk[-1]+d_fft[-1]))
    dfk=diff(edges)
    var0=(TI_ref*(0.75 + 5.6/uhub)*uhub)**2     # specified variance
    s0=var0-sum(G*dfk)                        # missing variance
    Sres=ones(size(dfk))*s0/size(dfk)        # spread missing variance over all frequencies
    # the discrete spectrum (see notebook2015, p. 23)
    S0=G*dfk+Sres
    return S0,fk,edges


def getDinstMatrix (P):
    n = P.shape[1]
    n=int(n)
    R=zeros((n,n))
    for ii in range(0,n):
        for jj in range(0,n):
            R[ii,jj]=linalg.norm(P[:,ii]-P[:,jj])
    return R

def CohFct (r,f,u,hubht):
    a=12                       #IEC 61400-3
    Lc=5.67*min(60,hubht)      # IEC 61400-3
    C=exp(-a*sqrt((f*r/u)**2+(0.12*r/Lc)**2))
    return C

def makeWindVeers (v0,P,Nf,dt_fft,plots,height,xi):
    Np=int(P.shape[1])
    hubht = height
    uhub  = v0
    TI_ref=0.16 # IEC normal turbulence model, turbine class A 
    # rhaghi, March 5, 2018: Changed to be in line with S class
    # turbulenceb!
    #rhaghi, March 5 2018: changed back to A class turb!
    #the spectrum
    [S_kai,fk,edges]= makeKaimalSpectrum (hubht, uhub, TI_ref,dt_fft,Nf)
    wk = 2*pi*fk
    T=ceil(1/fk[0])            # sampling window length
    t_fft=arange(0,T,dt_fft)
    R = getDinstMatrix (P)         # the distances between points
    # initializing variable
    S=zeros((int(R.shape[0]),int(R.shape[1]),Nf))
    H=zeros((int(R.shape[0]),int(R.shape[1]),Nf))
    Uk=1j*zeros((Np,Nf))
    u_fft_comp = 1j*zeros((t_fft.size,Np))
    u_fft = zeros((t_fft.size,Np))
    X=1j*zeros((Np,Np))
    phi=2*pi*xi
    
    for mm in range(0,Nf): ## the spectral matrix
        for kk in range(0,Np):
            for jj in range(0,Np):
                #print(mm,kk,jj)
                S[kk,jj,mm]=CohFct(R[kk,jj],fk[mm],uhub,hubht)*S_kai[mm]
        H[:,:,mm]=linalg.cholesky(abs(around(S[:,:,mm],decimals=10)))
        for kk in range(0,Np):
            X[kk,kk] =exp(1j*phi[kk,mm])
        #Fc = (H[:,:,mm]@X) @ ones((Np,1))
        Fc = matmul(matmul(H[:,:,mm],X),ones((Np,1)))
        Uk[:,mm] = Fc.reshape(1,Np)     ## get the Fourier coefficients wit random phases

## The iDFT with full random matrix (Veers' method)
#     fprintf('iDFT Veers\n')
    mm=0
    if plots ==1:
        pyplot.figure()
    
    for kk in range(0,Np):
        for tt in range(0,int(t_fft.size)):
            u_fft_comp[tt,kk] = sum((Uk[kk,:])/sqrt(2)*exp(1j*wk*t_fft[tt])+conj(Uk[kk,:]/sqrt(2)*exp(1j*wk*t_fft[tt])))
        u_fft[:,kk]=real(u_fft_comp[:,kk])+uhub;
        if plots ==1:
            pyplot.plot(t_fft,u_fft[:,mm])
            mm = mm+1
            pyplot.xlabel('Time [s]')
            pyplot.ylabel('Wind Speed [m/s]')
            pyplot.grid()
    #u_fft = u_fft.T
    return u_fft,t_fft,Uk,wk,edges,xi,S_kai,Fc,X

def makeWindVeers_630s (v0,P,Nf,dt_fft,plots,height,xi):
    Np=int(P.shape[1])
    hubht = height
    uhub  = v0
    TI_ref=0.16  # IEC normal turbulence model, turbine class A 
    # rhaghi, March 5, 2018: Changed to be in line with S class
    # turbulenceb!
    #rhaghi, March 5 2018: changed back to A class turb!
    #the spectrum
    [S_kai,fk,edges]= makeKaimalSpectrum_630s(hubht, uhub, TI_ref,dt_fft,Nf)
    wk = 2*pi*fk
    T=ceil(1/fk[0])            # sampling window length
    t_fft=arange(0,T,dt_fft)
    R = getDinstMatrix (P)         # the distances between points
    # initializing variable
    S=zeros((int(R.shape[0]),int(R.shape[1]),Nf))
    H=zeros((int(R.shape[0]),int(R.shape[1]),Nf))
    Uk=1j*zeros((Np,Nf))
    u_fft_comp = 1j*zeros((t_fft.size,Np))
    u_fft = zeros((t_fft.size,Np))
    X=1j*zeros((Np,Np))
    phi=2*pi*xi
    
    for mm in range(0,Nf): ## the spectral matrix
        for kk in range(0,Np):
            for jj in range(0,Np):
                #print(mm,kk,jj)
                S[kk,jj,mm]=CohFct(R[kk,jj],fk[mm],uhub,hubht)*S_kai[mm]
        H[:,:,mm]=linalg.cholesky(abs(around(S[:,:,mm],decimals=10)))
        for kk in range(0,Np):
            X[kk,kk] =exp(1j*phi[kk,mm])
        #Fc = (H[:,:,mm]@X) @ ones((Np,1))
        Fc = matmul(matmul(H[:,:,mm],X),ones((Np,1)))
        Uk[:,mm] = Fc.reshape(1,Np)     ## get the Fourier coefficients wit random phases

## The iDFT with full random matrix (Veers' method)
#     fprintf('iDFT Veers\n')
    mm=0
    if plots ==1:
        pyplot.figure()
    
    for kk in range(0,Np):
        for tt in range(0,int(t_fft.size)):
            u_fft_comp[tt,kk] = sum((Uk[kk,:])/sqrt(2)*exp(1j*wk*t_fft[tt])+conj(Uk[kk,:]/sqrt(2)*exp(1j*wk*t_fft[tt])))
        u_fft[:,kk]=real(u_fft_comp[:,kk])+uhub;
        if plots ==1:
            pyplot.plot(t_fft,u_fft[:,mm])
            mm = mm+1
            pyplot.xlabel('Time [s]')
            pyplot.ylabel('Wind Speed [m/s]')
            pyplot.grid()
    #u_fft = u_fft.T
    return u_fft,t_fft,Uk,wk,edges,xi,S_kai,Fc,X