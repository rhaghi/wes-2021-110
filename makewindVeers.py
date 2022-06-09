# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 17:39:59 2019

@author: radha
"""

## this function generates wind according to Veers model on the points provided in P
## based on severalPoints.m @rev 102
## 26 Feb 2018, Rad Haghi: This function upedate to have plot control and height as the input. Also, as we want to keep the randomness
## constant for all the trajectory, xi has taken out of the function to get realized once for wach realization on a trajectory

def makeWindVeers (v0,P,Nf,dt_fft,plots,height,xi):
#addpath('../../../Matlab_functionLib/')
Np=P[:,0].shape

## controls
#plots=plot;        # toggle plots on/ off

# fprintf('generating wind\n')


## parameters
hubht = height;
uhub  = v0;
TI_ref=0.16;    # IEC normal turbulence model, turbine class A 
                # rhaghi, March 5, 2018: Changed to be in line with S class
                # turbulence!
                #rhaghi, March 5 2018: changed back to A class turb!

## the spectrum
[S_kai,fk,edges]= makeKaimalSpectrum (hubht, uhub, TI_ref,dt_fft,Nf);

wk = 2*pi()*fk;
T=ceil(1/fk(1));            # sampling window length
t_fft=0:dt_fft:T;

# figure; hold off
# semilogx (fk,S0,'k-','LineWidth',2); hold on;


## generate matices for correlation
R = getDistMatrix (P);          # the distances between points

# initializing variables
S=zeros([size(R),Nf]);
H=zeros([size(R),Nf]);
Uk=zeros(Np,Nf);

u_fft  =zeros([length(t_fft),Np]);
X=zeros(Np);

if plots
    figure (2); leg2={}; hold off
end

# plot(edges,zeros(length(edges)),'o')
# plot(fk,zeros(length(fk)),'x')

# mar={':kx',':ks',':k*'};

# fprintf('calculating wind realization \n')
# the uncorrelated phases
#     rng(1);     #fix seed for one specific realization of all phases
rng('shuffle')      #shuffle to get one random realization of first phase
#xi=rand(Np,Nf);
phi=2*pi*xi;

for mm=1:Nf
    ## the spectral matrix
    for kk=1:Np
        for jj=1:Np
            S(kk,jj,mm)=CohFct(R(kk,jj),fk(mm),uhub,hubht)*S_kai(mm);    # the whole spectral matrix, Veers (2.1)
        end                                                     # homogenious wind, -> ss only a funcfion of fk, not of positon
    end
    
    ## solving for weighting matrix H (Veers1988, eq. 2.3)
    # see nb2015, p. 139
    for kk=1:Np
        for jj=1:kk-1
            H(kk,jj,mm)=(S(kk,jj,mm) - sum(H(kk,1:jj-1,mm).*H(jj,1:jj-1,mm))) / H(jj,jj,mm);
        end
        jj=kk;
        H(jj,jj,mm)=(S(jj,jj,mm) - sum(H(jj,1:jj-1,mm).^2))^(1/2);
    end
    
    
    ## get the Fourier coefficients wit random phases
    for kk=1:Np
        X(kk,kk)=exp(1i*phi(kk,mm));
    end
    
    Uk(:,mm)=H(:,:,mm)*X*ones(Np,1);
end

## The iDFT with full random matrix (Veers' method)
#     fprintf('iDFT Veers\n')
mm=1;
for kk=1:length(P(:,1))
    for tt=1:length(t_fft)
        u_fft(tt,kk)=sum((Uk(kk,:))/sqrt(2) .*exp( 1i*wk*t_fft(tt))+...
            conj((Uk(kk,:))/sqrt(2)).*exp(-1i*wk*t_fft(tt)));
    end
    u_fft(:,kk)=u_fft(:,kk)+uhub;
    
    if plots
        # visualize FFT data
        figure(2);
        plot(t_fft,u_fft(:,kk),'LineWidth',2); hold on;
        leg2=[leg2,sprintf('point (%.f, %.f) m,',P(kk,:))];
        #     drawnow
        mm=mm+1;
    end
end

u=u_fft'; 

 return u,t_fft,Uk,wk,edges,xi,S_kai
