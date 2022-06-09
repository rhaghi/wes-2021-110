# -*- coding: utf-8 -*-
"""
Created on Sat Aug  7 13:37:11 2021

@author: radha
"""
import chaospy as cp
import numpy as np
import matplotlib.pyplot as plt
import pickle
#from makewindVeers_main import *
#import seaborn as sns
import scipy.stats as st
import statistics as sts
import scipy.fftpack
import scipy.special
#from aP import *
#from HermitePoly import myHermite
import scipy.io as sio
from datetime import date
import time
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection # New import
#from LLFunc import *
import bottleneck as bn
import numpoly 
import fnmatch
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size = 14)
plt.rc('xtick', labelsize = 12)
plt.rc('ytick', labelsize = 12)

#%%
f = open('NREL5MW_trt_12mps_48000seeds.pckl','rb')
Trt_48K = pickle.load(f)
f.close()

Trt_48K_Mean = np.mean(Trt_48K,axis=0)
Trt_48K_Std = np.std(Trt_48K,axis=0)
Trt_48K_Q = np.zeros((3,Trt_48K_Mean.size))

Trt_48K_Q[0,:] = np.quantile(Trt_48K, 0.25,axis=0)
Trt_48K_Q[1,:] = np.quantile(Trt_48K, 0.5,axis=0)
Trt_48K_Q[2,:] = np.quantile(Trt_48K, 0.75,axis=0)


#%%
ITStepSecond = 30
ITStep = ITStepSecond*10
NoOfRlz = np.array([132,536,2002,6006]) #This is based on the polynomial order and the number of random variables
SimLength = np.array([27.3,10,10,10])*10+ITStep
SimLength_xaxis = np.array([27.3,10,10,10])*10
plot_indc=((0,0),(0,1),(1,0),(1,1))
PCEsimLength=np.array([27.3,6.8,1.8,1])
PolyOrder = np.array([2,3,4,5])
PCE_percent = np.array([25,50,75])

Trt_PCE_percentile = np.zeros((PolyOrder.size,600,PCE_percent.size))
dist = cp.Iid(cp.Uniform(0,1),10)
#%%
cc = fnmatch.filter(os.listdir(), "*_Thrust_*.np?")
cc.sort()
#%%
for i in range(0,SimLength.size):
    poly_load_thrust = cp.load(cc[i+1],allow_pickle=True)
    Trt_PCEApprox = poly_load_thrust['arr_0'].item()
    for j in range(ITStep,SimLength[i].astype(int)):
        for k in range(0,PCE_percent.size):
            PCE_poly = Trt_PCEApprox[str(0),str(j)]
            Q = cp.Perc(PCE_poly, PCE_percent[k] , dist,sample=48000)
            Trt_PCE_percentile[i,j,k] = Q
            print(i,j,k)

# tstr_file = time.strftime('%Y%m%d')


# f = open(tstr_file + '_Thrust_PCE_Percentile.pckl','wb')
# pickle.dump(Trt_PCE_percentile, f)
# f.close()
#%%

ITStepSecond = 30
ITStep = ITStepSecond*10
NoOfRlz = np.array([132,536,2002,6006]) #This is based on the polynomial order and the number of random variables
SimLength = np.array([27.3,10,10,10])*10+ITStep
SimLength_xaxis = np.array([27.3,10,10,10])*10
plot_indc=((0,0),(0,1),(1,0),(1,1))
PCEsimLength=np.array([27.3,6.8,1.8,1])
PolyOrder = np.array([2,3,4,5])
PCE_percent = np.array([25,50,75])

Trq_PCE_percentile = np.zeros((PolyOrder.size,600,PCE_percent.size))
dist = cp.Iid(cp.Uniform(0,1),10)

cc = fnmatch.filter(os.listdir(), "*_Torque_*.np?")
cc.sort()

for i in range(0,SimLength.size):
    poly_load_trq = cp.load(cc[i+1],allow_pickle=True)
    Trq_PCEApprox = poly_load_trq['arr_0'].item()
    for j in range(ITStep,SimLength[i].astype(int)):
        for k in range(0,PCE_percent.size):
            PCE_poly = Trq_PCEApprox[str(0),str(j)]
            Q = cp.Perc(PCE_poly, PCE_percent[k] , dist,sample=48000)
            Trq_PCE_percentile[i,j,k] = Q
            print(i,j,k)
            

# tstr_file = time.strftime('%Y%m%d')

# f = open(tstr_file + '_Torque_PCE_Percentile.pckl','wb')
# pickle.dump(Trq_PCE_percentile, f)
# f.close()
#%%

plt.close('all')    

visib=0.2
fig3, axes = plt.subplots(2,2,figsize=(11.69,8.27),sharey=True)


ITStepSecond = 30
ITStep = ITStepSecond*10
NoOfRlz = np.array([132,536,2002,6006]) #This is based on the polynomial order and the number of random variables
SimLength = np.array([27.3,10,10,10])*10+ITStep
SimLength_xaxis = np.array([27.3,10,10,10])*10
plot_indc=((0,0),(0,1),(1,0),(1,1))
PCEsimLength=np.array([27.3,6.8,1.8,1])
PolyOrder = np.array([2,3,4,5])
PCE_percent = np.array([25,50,75])

for i in range(0,NoOfRlz.size):
    x= np.arange(0,SimLength_xaxis[i].astype(int))/10
    axes[plot_indc[i]].clear()
    axes[plot_indc[i]].plot(x,Trt_PCE_percentile[i,ITStep:SimLength[i].astype(int),1]/1000,color='red',linewidth=3,label='PCE $Q_{2}$')
    axes[plot_indc[i]].plot(x,Trt_PCE_percentile[i,ITStep:SimLength[i].astype(int),0]/1000,'--',color='red',linewidth=3,label='PCE $Q_1$, $Q_3$')
    axes[plot_indc[i]].plot(x,Trt_PCE_percentile[i,ITStep:SimLength[i].astype(int),2]/1000,'--',color='red',linewidth=3)
    axes[plot_indc[i]].plot(x,Trt_48K_Q[1,ITStep:SimLength[i].astype(int)]/1000,color='blue',linewidth=3,label='48K $Q_2$')
    axes[plot_indc[i]].plot(x,Trt_48K_Q[0,ITStep:SimLength[i].astype(int)]/1000,'--',color='blue',linewidth=3,label='48K $Q_1$, $Q_3$')
    axes[plot_indc[i]].plot(x,Trt_48K_Q[2,ITStep:SimLength[i].astype(int)]/1000,'--',color='blue',linewidth=3)
    axes[plot_indc[i]].set_xlabel('Time [$s$]',fontsize=14)
    axes[plot_indc[i]].set_ylabel('Thrust [$kN$]',fontsize=14)
    axes[plot_indc[i]].grid()
    axes[plot_indc[i]].set_title('No of simulations: '+str(NoOfRlz[i].astype(int))+
        ', P='+str(PolyOrder[i].astype(int)), 
        horizontalalignment='center', fontsize=14)
    axes[plot_indc[i]].axvline(x=PCEsimLength[i],linewidth=2,ls=':',color='black')
    axes[plot_indc[i]].text((PCEsimLength[i]-0.1*PCEsimLength[i])/(x[-1]), 0.5, 'Cum Sim Lengh : 3600 Sec', rotation=90, verticalalignment='center'
        ,fontsize=14,transform= axes[plot_indc[i]].transAxes,)
    axes[plot_indc[i]].legend(fontsize = 12)
fig3.subplots_adjust(hspace=0.35)


tstr_fig = time.strftime('%Y%m%d')

fig3.savefig(tstr_fig+'Trt_PCE_SimsTimeSerQ1Q2Q3.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig3.savefig(tstr_fig+'Trt_PCE_SimsTimeSerQ1Q2Q3.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%
f = open('NREL5MW_trq_12mps_48000seeds.pckl','rb')
Trq_48K = pickle.load(f)
f.close()

Trq_48K_Mean = np.mean(Trq_48K,axis=0)
Trq_48K_Std = np.std(Trq_48K,axis=0)
Trq_48K_Q = np.zeros((3,Trq_48K_Mean.size))

Trq_48K_Q[0,:] = np.quantile(Trq_48K, 0.25,axis=0)
Trq_48K_Q[1,:] = np.quantile(Trq_48K, 0.5,axis=0)
Trq_48K_Q[2,:] = np.quantile(Trq_48K, 0.75,axis=0)


plt.close('all')    

visib=0.2
fig3, axes = plt.subplots(2,2,figsize=(11.69,8.27),sharey=True)


ITStepSecond = 30
ITStep = ITStepSecond*10
NoOfRlz = np.array([132,536,2002,6006]) #This is based on the polynomial order and the number of random variables
SimLength = np.array([27.3,10,10,10])*10+ITStep
SimLength_xaxis = np.array([27.3,10,10,10])*10
plot_indc=((0,0),(0,1),(1,0),(1,1))
PCEsimLength=np.array([27.3,6.8,1.8,1])
PolyOrder = np.array([2,3,4,5])
PCE_percent = np.array([25,50,75])

for i in range(0,NoOfRlz.size):
    x= np.arange(0,SimLength_xaxis[i].astype(int))/10
    axes[plot_indc[i]].clear()
    axes[plot_indc[i]].plot(x,Trq_PCE_percentile[i,ITStep:SimLength[i].astype(int),1]/1000,color='red',linewidth=3,label='PCE $Q_{2}$')
    axes[plot_indc[i]].plot(x,Trq_PCE_percentile[i,ITStep:SimLength[i].astype(int),0]/1000,'--',color='red',linewidth=3,label='PCE $Q_1$, $Q_3$')
    axes[plot_indc[i]].plot(x,Trq_PCE_percentile[i,ITStep:SimLength[i].astype(int),2]/1000,'--',color='red',linewidth=3)
    axes[plot_indc[i]].plot(x,Trq_48K_Q[1,ITStep:SimLength[i].astype(int)]/1000,color='blue',linewidth=3,label='48K $Q_2$')
    axes[plot_indc[i]].plot(x,Trq_48K_Q[0,ITStep:SimLength[i].astype(int)]/1000,'--',color='blue',linewidth=3,label='48K $Q_1$, $Q_3$')
    axes[plot_indc[i]].plot(x,Trq_48K_Q[2,ITStep:SimLength[i].astype(int)]/1000,'--',color='blue',linewidth=3)
    axes[plot_indc[i]].set_xlabel('Time [$s$]',fontsize=14)
    axes[plot_indc[i]].set_ylabel('Torque [$kNm$]',fontsize=14)
    axes[plot_indc[i]].grid()
    axes[plot_indc[i]].set_title('No of simulations: '+str(NoOfRlz[i].astype(int))+
        ', P='+str(PolyOrder[i].astype(int)), 
        horizontalalignment='center', fontsize=14)
    axes[plot_indc[i]].axvline(x=PCEsimLength[i],linewidth=2,ls=':',color='black')
    axes[plot_indc[i]].text((PCEsimLength[i]-0.1*PCEsimLength[i])/(x[-1]), 0.5, 'Cum Sim Lengh : 3600 Sec', rotation=90, verticalalignment='center'
        ,fontsize=14,transform= axes[plot_indc[i]].transAxes,)
    axes[plot_indc[i]].legend(fontsize = 12)
fig3.subplots_adjust(hspace=0.35)


tstr_fig = time.strftime('%Y%m%d')

fig3.savefig(tstr_fig+'Trq_PCE_SimsTimeSerQ1Q2Q3.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig3.savefig(tstr_fig+'Trq_PCE_SimsTimeSerQ1Q2Q3.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%


f1 = open('20210801_Thrust_MC_sims_P2_1M_rlz.pckl','rb')
f2 = open('20210801_Thrust_MC_sims_P3_1M_rlz.pckl','rb')
f3 = open('20210731_Thrust_MC_sims_P4_1M_rlz.pckl','rb')
f4 = open('20210726_Thrust_MC_sims_P5_1M_rlz.pckl','rb')
f5 = open('NREL5MW_trt_12mps_48000seeds.pckl','rb')

Trt_MC = np.zeros((4,1000000))

Trt_MC[0,:] = pickle.load(f1)
Trt_MC[1,:] = pickle.load(f2)
Trt_MC[2,:] = pickle.load(f3)
Trt_MC[3,:] = pickle.load(f4)
Trt_48K = pickle.load(f5)

f1.close()
f2.close()
f3.close()
f4.close()


mean = np.mean(Trt_MC)
standard_deviation = np.std(Trt_MC)
distance_from_mean = abs(Trt_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
Trt_MC_no_ol = Trt_MC[not_outlier]

#Trt_MC = Trt_MC/1000
#Trt_48K = Trt_48K/1000

Ntimes = Trt_48K.shape[1]

lb = np.min(Trt_MC_no_ol)
ub = np.max(Trt_MC_no_ol)

bins_o = np.linspace(lb,ub,30)

Trt_hist = np.zeros((bins_o.shape[0]-1,Ntimes))
Trt_MC_hist = np.zeros((bins_o.shape[0]-1,Trt_MC.shape[0]))
HL_dist= np.zeros((Trt_MC.shape[0],Ntimes))
Trt_hist = np.zeros(bins_o.shape[0]-1)
Trt_MC_hist = np.zeros((bins_o.shape[0]-1,Trt_MC.shape[0]))



for i in range(0,Trt_MC.shape[0]):
    #lb = np.min(Trt_MC[i,:])
    #ub = np.max(Trt_MC[i,:])
    bins_i = np.linspace(lb,ub,30)
    Trt_MC_hist[:,i],_ = np.histogram(Trt_MC[i,:],bins=bins_i,density='True')
    Q = np.sqrt(Trt_MC_hist*np.diff(bins_i)[0])
    #Q = Trt_MC_hist
    for j in range(0,Ntimes):    
        Trt_hist,blabla = np.histogram(Trt_48K[:,j],bins=bins_i,density='True')
        P = np.sqrt(Trt_hist*np.diff(bins_i)[0])
        #P = Trt_hist
        HH = (P - Q[:,i])**2
        HL_dist[i,j] = (1/np.sqrt(2))*(np.sqrt(np.sum(HH)))
        #HL_dist[i,j] = st.entropy(P, Q[:,i])
        print(i,j)


timestr = time.strftime('%Y%m%d')  

f=open(timestr+'_HellingerDist_Trt_vs_MC.pckl','wb')
pickle.dump(HL_dist,f)
f.close()    

timestr = time.strftime('%Y%m%d')  
plt.close('all')
t_fft = np.linspace(0,600,Ntimes)

plt.close('all')    

visib=0.3
fig, axes = plt.subplots(1,1,figsize=(11.69,8.27),sharey=True)
PolyOrder = np.array([2,3,4,5])


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size = 18)
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)

for i in range(0,Trt_MC.shape[0]):
    axes.plot(t_fft,HL_dist[i,:]*100,'o',label = 'P = '+str(PolyOrder[i].astype(int)))
    axes.set_xlabel('Time [s]')
    axes.set_ylabel('Thrust Hellinger distance [\%]')
    #axes.set_yticks(np.arange(0,3.25, step=0.25))

axes.xaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.yaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.xaxis.grid(which='minor',lw=0.75,alpha=0.5)
axes.yaxis.grid(which='minor',lw=0.75,alpha=0.5)
#axes.grid(True, which='both')
axes.minorticks_on()
axes.legend()

fig.savefig(timestr+'_HL_Dis_time_Trt_vs_MCs.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr+'_HL_Dis_time_Trt_vs_MCs.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)




f1 = open('20210801_Torque_MC_sims_P2_1M_rlz.pckl','rb')
f2 = open('20210801_Torque_MC_sims_P3_1M_rlz.pckl','rb')
f3 = open('20210731_Torque_MC_sims_P4_1M_rlz.pckl','rb')
f4 = open('20210726_Torque_MC_sims_P5_1M_rlz.pckl','rb')
f5 = open('NREL5MW_trq_12mps_48000seeds.pckl','rb')

Trq_MC = np.zeros((4,1000000))

Trq_MC[0,:] = pickle.load(f1)
Trq_MC[1,:] = pickle.load(f2)
Trq_MC[2,:] = pickle.load(f3)
Trq_MC[3,:] = pickle.load(f4)
Trq_48K = pickle.load(f5)

f1.close()
f2.close()
f3.close()
f4.close()


Ntimes = Trq_48K.shape[1]

lb = np.min(Trq_48K)
ub = np.max(Trq_48K)

bins_o = np.linspace(lb,ub,30)

Trq_hist = np.zeros((bins_o.shape[0]-1,Ntimes))
Trq_MC_hist = np.zeros((bins_o.shape[0]-1,Trq_MC.shape[0]))
HL_dist= np.zeros((Trq_MC.shape[0],Ntimes))
Trq_hist = np.zeros(bins_o.shape[0]-1)
Trq_MC_hist = np.zeros((bins_o.shape[0]-1,Trq_MC.shape[0]))

for i in range(0,Trq_MC.shape[0]):
    #lb = np.min(Trq_MC[i,:])
    #ub = np.max(Trq_MC[i,:])
    bins_i = np.linspace(lb,ub,30)
    Trq_MC_hist[:,i],_ = np.histogram(Trq_MC[i,:],bins=bins_i,density='True')
    Q = np.sqrt(Trq_MC_hist*np.diff(bins_i)[0])
    #Q = Trq_MC_hist
    for j in range(0,Ntimes):    
        Trq_hist,blabla = np.histogram(Trq_48K[:,j],bins=bins_i,density='True')
        P = np.sqrt(Trq_hist*np.diff(bins_i)[0])
        #P = Trq_hist
        HH = (P - Q[:,i])**2
        HL_dist[i,j] = (1/np.sqrt(2))*(np.sqrt(np.sum(HH)))
        #HL_dist[i,j] = st.entropy(P, Q[:,i])
        print(i,j)

timestr = time.strftime('%Y%m%d')  

f=open(timestr+'_HellingerDist_Trq_vs_MC.pckl','wb')
pickle.dump(HL_dist,f)
f.close()    


timestr = time.strftime('%Y%m%d')  
plt.close('all')
t_fft = np.linspace(0,600,Ntimes)

plt.close('all')    

visib=0.3
fig, axes = plt.subplots(1,1,figsize=(11.69,8.27),sharey=True)
PolyOrder = np.array([2,3,4,5])


plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size = 18)
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)

for i in range(0,Trq_MC.shape[0]):
    axes.plot(t_fft,HL_dist[i,:]*100,'o',label = 'P = '+str(PolyOrder[i].astype(int)))
    axes.set_xlabel('Time [s]')
    axes.set_ylabel('Torque Hellinger distance [\%]')
    #axes.set_yticks(np.arange(0,3.25, step=0.25))

axes.xaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.yaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.xaxis.grid(which='minor',lw=0.75,alpha=0.5)
axes.yaxis.grid(which='minor',lw=0.75,alpha=0.5)
#axes.grid(True, which='both')
axes.minorticks_on()
axes.legend()

fig.savefig(timestr+'_HL_Dis_time_Trq_vs_MCs.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr+'_HL_Dis_time_Trq_vs_MCs.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)
#%%
f1 = open('20210801_Thrust_MC_sims_P2_1M_rlz.pckl','rb')
f2 = open('20210801_Thrust_MC_sims_P3_1M_rlz.pckl','rb')
f3 = open('20210731_Thrust_MC_sims_P4_1M_rlz.pckl','rb')
f4 = open('20210726_Thrust_MC_sims_P5_1M_rlz.pckl','rb')
f5 = open('NREL5MW_trt_12mps_48000seeds.pckl','rb')

Trt_MC = np.zeros((4,1000000))

Trt_MC[0,:] = pickle.load(f1)
Trt_MC[1,:] = pickle.load(f2)
Trt_MC[2,:] = pickle.load(f3)
Trt_MC[3,:] = pickle.load(f4)
Trt_48K = pickle.load(f5)

f1.close()
f2.close()
f3.close()
f4.close()



Trt_MC = Trt_MC/1000
Trt_48K = Trt_48K/1000




mean = np.mean(Trt_MC)
standard_deviation = np.std(Trt_MC)
distance_from_mean = abs(Trt_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
Trt_MC_no_ol = Trt_MC[not_outlier]

#%%
plt.close('all')

plt.rc('font', size = 18)
plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)

visib=0.3
fig3, axes = plt.subplots(2,2,figsize=(18,12.73),sharey=True)
timestr = time.strftime('%Y%m%d')  



k=0

plot_indc = ()

for i in range(0,6):
    for j in range(0,2):
        k=k+1
        plot_indc = plot_indc+((i,j),)

ITStepSecond = 30
ITStep = ITStepSecond*10
NoOfRlz = np.array([132,536,2002,6006]) #This is based on the polynomial order and the number of random variables
SimLength = np.array([27.3,10,10,10])*10+ITStep
SimLength_xaxis = np.array([27.3,10,10,10])*10
PCEsimLength=np.array([27.3,6.8,1.8,1])
PolyOrder = np.array([2,3,4,5])
elt = np.array([9.942,40.682,144.025,445.721])
l=0;
#PCEIndex[i].astype(int)

Trt48KIndex = np.random.randint(60,5999)

lb = np.min(Trt_MC_no_ol)
ub = np.max(Trt_MC_no_ol)

bins = np.linspace(lb,ub,40)


for i in range(0,PolyOrder.size):
    hist_Trt_MC,bin_Trt_MC = np.histogram(Trt_MC[i,:],bins = bins,density=True)
    hist_Trt_48K,bin_Trt_48K = np.histogram(Trt_48K[:,Trt48KIndex],bin_Trt_MC,density=True)
    axes[plot_indc[l]].bar(bin_Trt_48K[:-1],hist_Trt_48K* np.diff(bin_Trt_48K)*100,width =np.diff(bin_Trt_48K)[0],
        edgecolor='black',label = '48K sims',alpha = 0.7)
    axes[plot_indc[l]].bar(bin_Trt_MC[:-1],hist_Trt_MC* np.diff(bin_Trt_MC)*100,width =np.diff(bin_Trt_MC)[0],
        edgecolor='black',label='SM, 1e6 MCs, P='+str(PolyOrder[i].astype(int)),alpha = 0.7)
    axes[plot_indc[l]].legend()
        #axes[plot_indc[l]].set_xlabel('Bound circulation $\Gamma$ [$m^2/s$]')
    axes[plot_indc[l]].set_ylabel('Frequency of Occurances [\%]',fontsize=20)
    axes[plot_indc[l]].grid()
    #print(PCEIndex[j,i])
    l=l+1
        
axes[plot_indc[l-1]].set_xlabel('Thrust [kN]',fontsize=20)
axes[plot_indc[l-2]].set_xlabel('Thrust [kN]',fontsize=20)



fig3.savefig(timestr+'_MCPCEvsSim_Hist_Trt.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig3.savefig(timestr+'_MCPCEvsSim_Hist_Trt.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%

f1 = open('20210801_Torque_MC_sims_P2_1M_rlz.pckl','rb')
f2 = open('20210801_Torque_MC_sims_P3_1M_rlz.pckl','rb')
f3 = open('20210731_Torque_MC_sims_P4_1M_rlz.pckl','rb')
f4 = open('20210726_Torque_MC_sims_P5_1M_rlz.pckl','rb')
f5 = open('NREL5MW_trq_12mps_48000seeds.pckl','rb')

Trq_MC = np.zeros((4,1000000))

Trq_MC[0,:] = pickle.load(f1)
Trq_MC[1,:] = pickle.load(f2)
Trq_MC[2,:] = pickle.load(f3)
Trq_MC[3,:] = pickle.load(f4)
Trq_48K = pickle.load(f5)

f1.close()
f2.close()
f3.close()
f4.close()

Trq_MC = Trq_MC/1000
Trq_48K = Trq_48K/1000




mean = np.mean(Trq_MC)
standard_deviation = np.std(Trq_MC)
distance_from_mean = abs(Trq_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
Trq_MC_no_ol = Trq_MC[not_outlier]

#%%
plt.close('all')

plt.rc('font', size = 18)
plt.rc('xtick', labelsize = 18)
plt.rc('ytick', labelsize = 18)

visib=0.3
fig3, axes = plt.subplots(2,2,figsize=(18,12.73),sharey=True)
timestr = time.strftime('%Y%m%d')  



k=0

plot_indc = ()

for i in range(0,6):
    for j in range(0,2):
        k=k+1
        plot_indc = plot_indc+((i,j),)

ITStepSecond = 30
ITStep = ITStepSecond*10
NoOfRlz = np.array([132,536,2002,6006]) #This is based on the polynomial order and the number of random variables
SimLength = np.array([27.3,10,10,10])*10+ITStep
SimLength_xaxis = np.array([27.3,10,10,10])*10
PCEsimLength=np.array([27.3,6.8,1.8,1])
PolyOrder = np.array([2,3,4,5])
elt = np.array([9.942,40.682,144.025,445.721])
l=0;
#PCEIndex[i].astype(int)

Trq48KIndex = np.random.randint(60,5999)

lb = np.min(Trq_MC_no_ol)
ub = np.max(Trq_MC_no_ol)

bins = np.linspace(lb,ub,40)


for i in range(0,PolyOrder.size):
    hist_Trq_MC,bin_Trq_MC = np.histogram(Trq_MC[i,:],bins = bins,density=True)
    hist_Trq_48K,bin_Trq_48K = np.histogram(Trq_48K[:,Trq48KIndex],bin_Trq_MC,density=True)
    axes[plot_indc[l]].bar(bin_Trq_48K[:-1],hist_Trq_48K* np.diff(bin_Trq_48K)*100,width =np.diff(bin_Trq_48K)[0],
        edgecolor='black',label = '48K sims',alpha = 0.7)
    axes[plot_indc[l]].bar(bin_Trq_MC[:-1],hist_Trq_MC* np.diff(bin_Trq_MC)*100,width =np.diff(bin_Trq_MC)[0],
        edgecolor='black',label='SM, 1e6 MCs, P='+str(PolyOrder[i].astype(int)),alpha = 0.7)
    axes[plot_indc[l]].legend()
        #axes[plot_indc[l]].set_xlabel('Bound circulation $\Gamma$ [$m^2/s$]')
    axes[plot_indc[l]].set_ylabel('Frequency of Occurances [\%]',fontsize=20)
    axes[plot_indc[l]].grid()
    #print(PCEIndex[j,i])
    l=l+1
        
axes[plot_indc[l-1]].set_xlabel('Torque [kNm]',fontsize=20)
axes[plot_indc[l-2]].set_xlabel('Torque [kNm]',fontsize=20)



fig3.savefig(timestr+'_MCPCEvsSim_Hist_Trq.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig3.savefig(timestr+'_MCPCEvsSim_Hist_Trq.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%

f = open('NREL5MW_trt_12mps_48000seeds.pckl','rb')
Trt_48K =  pickle.load(f)
f.close()

a = np.arange(0,48000)
bb = np.random.permutation(a)
cc = bb.reshape((8000,6))
Trt_48K_g_stats = np.zeros((8000,5))
Trt_group = np.zeros((6,6000))

Trt_max_prct_48K = np.zeros(4)
Trt_max_prct_48K[0] = np.max(Trt_48K)
Trt_max_prct_48K[1] = np.percentile(Trt_48K,99)
Trt_max_prct_48K[2] = np.percentile(Trt_48K,95)
Trt_max_prct_48K[3] = np.percentile(Trt_48K,90)


for i in range(0,cc.shape[0]):
    Trt_group = Trt_48K[cc[i,:],:]
    Trt_48K_g_stats[i,0] = np.max(Trt_group)
    Trt_48K_g_stats[i,1] = np.mean(np.max(Trt_group,axis=1))
    Trt_48K_g_stats[i,2] = np.percentile(Trt_group,99)
    Trt_48K_g_stats[i,3] = np.percentile(Trt_group,95)
    Trt_48K_g_stats[i,4] = np.percentile(Trt_group,90)

#%%

MC_sims = ['10','100','1k','10k','48k','100k','1M','10M','100M','288M','500M']
time_str = '20210731'
file_name_temp = time_str+'_Thrust_MC_sims_P4_'

nbins = 30
Trt_max_prct = np.zeros((len(MC_sims),4))
Trt_bin_edge = np.zeros((nbins+1,len(MC_sims)))
Trt_bin_hst = np.zeros((nbins,len(MC_sims)))
i=0


for x in MC_sims:
    print(x)
    f = open(file_name_temp + x + '_rlz.pckl','rb')
    Trt = pickle.load(f)
    f.close()
    Trt_max_prct[i,0] = np.max(Trt)
    Trt_max_prct[i,1] = np.percentile(Trt,99)
    Trt_max_prct[i,2] = np.percentile(Trt,95)
    Trt_max_prct[i,3] = np.percentile(Trt,90)
    Trt_bin_hst[:,i] , Trt_bin_edge[:,i] =np.histogram(Trt,bins=nbins,density=True)
    i = i+1
    del Trt
#%%
plt.close('all')    
timestr_fig = time.strftime('%Y%m%d')  
plt.rc('font', size = 22)
plt.rc('xtick', labelsize = 22)
plt.rc('ytick', labelsize = 22)
visib=0.3
fig, axes = plt.subplots(1,1,figsize=(16.96,12),sharey=True)



MC_sims_num = np.array([10,100,1000,1e4,4.8e4,1e5,1e6,1e7,1e8,288e6,500e6])

cl = ['tab:blue','tab:orange','tab:green','tab:red']
lbls = ['Max','$P_{99}$','$P_{95}$','$P_{90}$']
lbls_48K = ['48K Max','48K $P_{99}$','48K $P_{95}$','48K $P_{90}$']
lbls_aSim = ['6 sims Max mean','6 sims $P_{99}$','6 sims $P_{95}$','6 sims $P_{90}$']





for i in range(0,Trt_max_prct.shape[1]):
     axes.semilogx(MC_sims_num,Trt_max_prct[:,i]/1000,color = cl[i],lw = 4,label = lbls[i])

for i in range(0,Trt_max_prct.shape[1]):
    axes.hlines(Trt_max_prct_48K[i]/1000,np.min(MC_sims_num),np.max(MC_sims_num),lw = 4, linestyles = (0, (1, 1)) ,color = cl[i],label = lbls_48K[i] )

for i in range(0,Trt_max_prct.shape[1]):
     axes.plot(np.nan,'o',color = cl[i],label = lbls_aSim [i])   

axes.legend(ncol = 1,loc='right',bbox_to_anchor=(1.27,0.5))



axes.set_xlim((6,7.5e8))
#axes.set_ylim((300,1200))
axes.xaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.yaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.xaxis.grid(which='minor',lw=0.75,alpha=0.5)
axes.yaxis.grid(which='minor',lw=0.75,alpha=0.5)
#axes.grid(True, which='both')
axes.minorticks_on()
box = axes.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

axes.set_xlabel('No. MC simulations [-]')
axes.set_ylabel('Thrust [kN]')

timestr_fig = time.strftime('%Y%m%d')  

axes2 = axes.twiny()

lbls = ['Max','$P_{99}$','$P_{95}$','$P_{90}$']

cl = ['tab:blue','tab:orange','tab:green','tab:red']


for i in range(1,Trt_48K_g_stats.shape[1]):
    axes2.plot(Trt_48K_g_stats[:,i]/1000,'o',color = cl[i-1],label = lbls[i-1],alpha=0.03)
    
axes2.set_xlabel('Aerodynamic simulations group no[-]')

fig.tight_layout()

fig.savefig(timestr_fig+'_TRT_ExtrmeLoads_AsimG_MCvsRef_P5.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr_fig+'_TRT_ExtrmeLoads_AsimG_MCvsRef_P5.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)
#%%


f = open('NREL5MW_trq_12mps_48000seeds.pckl','rb')
Trq_48K =  pickle.load(f)
f.close()

a = np.arange(0,48000)
bb = np.random.permutation(a)
cc = bb.reshape((8000,6))
Trq_48K_g_stats = np.zeros((8000,5))
Trq_group = np.zeros((6,6000))

Trq_max_prct_48K = np.zeros(4)
Trq_max_prct_48K[0] = np.max(Trq_48K)
Trq_max_prct_48K[1] = np.percentile(Trq_48K,99)
Trq_max_prct_48K[2] = np.percentile(Trq_48K,95)
Trq_max_prct_48K[3] = np.percentile(Trq_48K,90)


for i in range(0,cc.shape[0]):
    Trq_group = Trq_48K[cc[i,:],:]
    Trq_48K_g_stats[i,0] = np.max(Trq_group)
    Trq_48K_g_stats[i,1] = np.mean(np.max(Trq_group,axis=1))
    Trq_48K_g_stats[i,2] = np.percentile(Trq_group,99)
    Trq_48K_g_stats[i,3] = np.percentile(Trq_group,95)
    Trq_48K_g_stats[i,4] = np.percentile(Trq_group,90)

#%%

MC_sims = ['10','100','1k','10k','48k','100k','1M','10M','100M','288M','500M']
time_str = '20210731'
file_name_temp = time_str+'_Torque_MC_sims_P4_'

nbins = 30
Trq_max_prct = np.zeros((len(MC_sims),4))
Trq_bin_edge = np.zeros((nbins+1,len(MC_sims)))
Trq_bin_hst = np.zeros((nbins,len(MC_sims)))
i=0


for x in MC_sims:
    print(x)
    f = open(file_name_temp + x + '_rlz.pckl','rb')
    Trq = pickle.load(f)
    f.close()
    Trq_max_prct[i,0] = np.max(Trq)
    Trq_max_prct[i,1] = np.percentile(Trq,99)
    Trq_max_prct[i,2] = np.percentile(Trq,95)
    Trq_max_prct[i,3] = np.percentile(Trq,90)
    Trq_bin_hst[:,i] , Trq_bin_edge[:,i] =np.histogram(Trq,bins=nbins,density=True)
    i = i+1
    del Trq
#%%
plt.close('all')    
timestr_fig = time.strftime('%Y%m%d')  
plt.rc('font', size = 22)
plt.rc('xtick', labelsize = 22)
plt.rc('ytick', labelsize = 22)
visib=0.3
fig, axes = plt.subplots(1,1,figsize=(16.96,12),sharey=True)



MC_sims_num = np.array([10,100,1000,1e4,4.8e4,1e5,1e6,1e7,1e8,288e6,500e6])

cl = ['tab:blue','tab:orange','tab:green','tab:red']
lbls = ['Max','$P_{99}$','$P_{95}$','$P_{90}$']
lbls_48K = ['48K Max','48K $P_{99}$','48K $P_{95}$','48K $P_{90}$']
lbls_aSim = ['6 sims Max mean','6 sims $P_{99}$','6 sims $P_{95}$','6 sims $P_{90}$']





for i in range(0,Trq_max_prct.shape[1]):
     axes.semilogx(MC_sims_num,Trq_max_prct[:,i]/1000,color = cl[i],lw = 4,label = lbls[i])

for i in range(0,Trq_max_prct.shape[1]):
    axes.hlines(Trq_max_prct_48K[i]/1000,np.min(MC_sims_num),np.max(MC_sims_num),lw = 4, linestyles = (0, (1, 1)) ,color = cl[i],label = lbls_48K[i] )

for i in range(0,Trq_max_prct.shape[1]):
     axes.plot(np.nan,'o',color = cl[i],label = lbls_aSim [i])   

axes.legend(ncol = 1,loc='right',bbox_to_anchor=(1.27,0.5))



axes.set_xlim((6,7.5e8))
#axes.set_ylim((300,1200))
axes.xaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.yaxis.grid(which='major',lw=1.5,alpha=0.5)
axes.xaxis.grid(which='minor',lw=0.75,alpha=0.5)
axes.yaxis.grid(which='minor',lw=0.75,alpha=0.5)
#axes.grid(True, which='both')
axes.minorticks_on()
box = axes.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

axes.set_xlabel('No. MC simulations [-]')
axes.set_ylabel('Torque [kNm]')

timestr_fig = time.strftime('%Y%m%d')  

axes2 = axes.twiny()

lbls = ['Max','$P_{99}$','$P_{95}$','$P_{90}$']

cl = ['tab:blue','tab:orange','tab:green','tab:red']


for i in range(1,Trq_48K_g_stats.shape[1]):
    axes2.plot(Trq_48K_g_stats[:,i]/1000,'o',color = cl[i-1],label = lbls[i-1],alpha=0.03)
    
axes2.set_xlabel('Aerodynamic simulations group no[-]')

fig.tight_layout()

fig.savefig(timestr_fig+'_Trq_ExtrmeLoads_AsimG_MCvsRef_P5.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr_fig+'_Trq_ExtrmeLoads_AsimG_MCvsRef_P5.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%
f=open('20220424_HellingerDist_Trt_vs_MC.pckl','rb')
trt_HL_MC_48k = pickle.load(f)
f.close()

f=open('20220424_HellingerDist_Trq_vs_MC.pckl','rb')
trq_HL_MC_48k = pickle.load(f)
f.close()
#%%

f = open('NREL5MW_trt_12mps_48000seeds.pckl','rb')
Trt_48K =  pickle.load(f)
f.close()

f=open('20210811_Thrust_PCE_Percentile.pckl','rb')
trt_pce_pct = pickle.load(f)
f.close()

Trt_48K_Mean = np.mean(Trt_48K,axis=0)
Trt_48K_Std = np.std(Trt_48K,axis=0)


Trt_48K_Q = np.zeros((3,Trt_48K_Mean.size))

Trt_48K_Q[0,:] = np.quantile(Trt_48K, 0.25,axis=0)
Trt_48K_Q[1,:] = np.quantile(Trt_48K, 0.5,axis=0)
Trt_48K_Q[2,:] = np.quantile(Trt_48K, 0.75,axis=0)
#%%
trt_q_err = np.zeros((2,3))
k=0
for i in range(2,4):
    trt_q_err[k,:] = np.mean(100*((Trt_48K_Q[:,300:400]-trt_pce_pct[i,300:400,:].T)/Trt_48K_Q[:,300:400]),axis=1)
    k=k+1


#%%

f = open('NREL5MW_trq_12mps_48000seeds.pckl','rb')
Trq_48K =  pickle.load(f)
f.close()

f=open('20210811_Torque_PCE_Percentile.pckl','rb')
trq_pce_pct = pickle.load(f)
f.close()


Trq_48K_Mean = np.mean(Trq_48K,axis=0)
Trq_48K_Std = np.std(Trq_48K,axis=0)


Trq_48K_Q = np.zeros((3,Trq_48K_Mean.size))

Trq_48K_Q[0,:] = np.quantile(Trq_48K, 0.25,axis=0)
Trq_48K_Q[1,:] = np.quantile(Trq_48K, 0.5,axis=0)
Trq_48K_Q[2,:] = np.quantile(Trq_48K, 0.75,axis=0)

#%%
trq_q_err = np.zeros((2,3))
k=0
for i in range(2,4):
    trq_q_err[k,:] = np.mean(100*((Trq_48K_Q[:,300:400]-trq_pce_pct[i,300:400,:].T)/Trq_48K_Q[:,300:400]),axis=1)
    k=k+1



#%%

f = open('NREL5MW_trt_12mps_48000seeds.pckl','rb')
Trt_48K =  pickle.load(f)
f.close()

Ntimes = Trt_48K.shape[1]

lb = np.min(Trt_48K)
ub = np.max(Trt_48K)

bins = np.arange(lb,ub,(ub-lb)/20)

Trt_hist = np.zeros((bins.shape[0]-1,Ntimes))
HL_dist= np.zeros((Ntimes,Ntimes))
    
for i in range(0,Ntimes):
    Trt_hist[:,i],blabla = np.histogram(Trt_48K[:,i],bins=bins,density='True')
    print(i)

P = np.sqrt(Trt_hist*np.diff(bins)[0])
#Q = np.sqrt(Gamma_un_hist)

for i in range(0,Ntimes):
    for j in range(0,Ntimes):
        HH = (P[:,i] - P[:,j])**2
        HL_dist[i,j] = (1/np.sqrt(2))*(np.sqrt(np.sum(HH,axis=0)))
        print(i,j)
        
timestr = time.strftime('%Y%m%d')  

#%%

f = open('NREL5MW_trt_12mps_48000seeds.pckl','rb')
Trt_48K =  pickle.load(f)
f.close()

Ntimes = Trt_48K.shape[1]


f=open('20220424_HellingerDist_Trt.pckl','rb')
HL_dist = pickle.load(f)
f.close()

timestr = time.strftime('%Y%m%d')  
plt.close('all')
t_fft = np.linspace(0,600,Ntimes)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size = 18)
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)



fig, ax = plt.subplots(1, 1,figsize=(11.69,8.27))

ax.plot(t_fft[25:],np.max(HL_dist[25:,25:],axis=0)*100,'o',alpha=0.2,label='Max. Hellinger Distance')
ub = np.max(np.max(HL_dist[25:,25:],axis=0)*100)
lb = np.min(np.max(HL_dist[25:,25:],axis=0)*100)
ax.grid()
ax.hlines(ub,0,600,lw=3,linestyles='dashed',color='C1',label='Upper bound= '+str(np.round(ub,2))+'\%')
ax.hlines(lb,0,600,lw=3,linestyles='dashed',color='C2',label='Lower bound= '+str(np.round(lb,2))+'\%')
box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(ncol = 3,loc='upper center', bbox_to_anchor=(0.5, 1.1))
#ax.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
ax.set_xlabel('Time [s]')
ax.set_ylabel('Max. Hellinger distance for Thrust [\%]')

fig.savefig(timestr+'_HL_Dis_time_Trt.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr+'_HL_Dis_time_Trt.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%

f = open('NREL5MW_trq_12mps_48000seeds.pckl','rb')
Trq_48K =  pickle.load(f)
f.close()

Ntimes = Trq_48K.shape[1]

lb = np.min(Trq_48K)
ub = np.max(Trq_48K)

bins = np.arange(lb,ub,(ub-lb)/20)

Trq_hist = np.zeros((bins.shape[0]-1,Ntimes))
HL_dist= np.zeros((Ntimes,Ntimes))
    
for i in range(0,Ntimes):
    Trq_hist[:,i],blabla = np.histogram(Trq_48K[:,i],bins=bins,density='True')
    print(i)

P = np.sqrt(Trq_hist*np.diff(bins)[0])
#Q = np.sqrt(Gamma_un_hist)

for i in range(0,Ntimes):
    for j in range(0,Ntimes):
        HH = (P[:,i] - P[:,j])**2
        HL_dist[i,j] = (1/np.sqrt(2))*(np.sqrt(np.sum(HH,axis=0)))
        print(i,j)
        
timestr = time.strftime('%Y%m%d')  

#%%
f = open('NREL5MW_trq_12mps_48000seeds.pckl','rb')
Trq_48K =  pickle.load(f)
f.close()

Ntimes = Trq_48K.shape[1]

f=open('20220424_HellingerDist_Trq.pckl','rb')
HL_dist = pickle.load(f)
f.close()
#%%
timestr = time.strftime('%Y%m%d')  
plt.close('all')
t_fft = np.linspace(0,600,Ntimes)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size = 18)
plt.rc('xtick', labelsize = 16)
plt.rc('ytick', labelsize = 16)



fig, ax = plt.subplots(1, 1,figsize=(11.69,8.27))

ax.plot(t_fft[25:],np.max(HL_dist[25:,25:],axis=0)*100,'o',alpha=0.2,label='Max. Hellinger Distance')
ub = np.max(np.max(HL_dist[25:,25:],axis=0)*100)
lb = np.min(np.max(HL_dist[25:,25:],axis=0)*100)
ax.grid()
ax.hlines(ub,0,600,lw=3,linestyles='dashed',color='C1',label='Upper bound= '+str(np.round(ub,2))+'\%')
ax.hlines(lb,0,600,lw=3,linestyles='dashed',color='C2',label='Lower bound= '+str(np.round(lb,2))+'\%')
box = ax.get_position()
#ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(ncol = 3,loc='upper center', bbox_to_anchor=(0.5, 1.1))
#ax.legend(loc='best', bbox_to_anchor=(0.5, 0., 0.5, 0.5))
ax.set_xlabel('Time [s]')
ax.set_ylabel('Max. Hellinger distance for Torque[\%]')

fig.savefig(timestr+'_HL_Dis_time_Trq.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr+'_HL_Dis_time_Trq.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%