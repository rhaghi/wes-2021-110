# -*- coding: utf-8 -*-
"""
Created on Thu Jan  3 10:05:44 2019

@author: radha
"""

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