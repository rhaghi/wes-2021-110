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
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # New import
#from LLFunc import *
import bottleneck as bn
import numpoly
import fnmatch
import os

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)
plt.rc('xtick', labelsize=12)
plt.rc('ytick', labelsize=12)

# %%


f2 = open('20220308_Thrust_MC_sims_P2_1M_rlz.pckl', 'rb')
f3 = open('20220308_Thrust_MC_sims_P3_1M_rlz.pckl', 'rb')
f4 = open('20220324_Thrust_MC_sims_P4_1M_rlz.pckl', 'rb')

f5 = open('NREL5MW_trt_12mps_48000seeds.pckl', 'rb')

Trt_MC = np.zeros((4, 1000000))

Trt_MC[0, :] = pickle.load(f2)
Trt_MC[1, :] = pickle.load(f3)
Trt_MC[2, :] = pickle.load(f4)

Trt_48K = pickle.load(f5)


f2.close()
f3.close()
f5.close()


Trt_MC = Trt_MC/1000
Trt_48K = Trt_48K/1000


mean = np.mean(Trt_MC)
standard_deviation = np.std(Trt_MC)
distance_from_mean = abs(Trt_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
Trt_MC_no_ol = Trt_MC[not_outlier]


plt.close('all')

plt.rc('font', size=18)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

visib = 0.3
fig3, axes = plt.subplots(1, 3, figsize=(27, 6.35), sharey=True)
timestr = time.strftime('%Y%m%d')


k = 0

plot_indc = ()

for i in range(0, 1):
    for j in range(0, 2):
        k = k+1
        plot_indc = plot_indc+((i, j),)

ITStepSecond = 30
ITStep = ITStepSecond*10
# This is based on the polynomial order and the number of random variables
NoOfRlz = np.array([132, 536, 2002, 6006])
SimLength = np.array([27.3, 10, 10, 10])*10+ITStep
SimLength_xaxis = np.array([27.3, 10, 10, 10])*10
PCEsimLength = np.array([27.3, 6.8, 1.8, 1])
PolyOrder = np.array([2, 3, 4])
elt = np.array([9.942, 40.682, 144.025, 445.721])
l = 0
# PCEIndex[i].astype(int)

Trt48KIndex = np.random.randint(60, 5999)

lb = np.min(Trt_MC_no_ol)
ub = np.max(Trt_MC_no_ol)

bins = np.linspace(lb, ub, 40)


for i in range(0, PolyOrder.size):
    hist_Trt_MC, bin_Trt_MC = np.histogram(
        Trt_MC[i, :], bins=bins, density=True)
    hist_Trt_48K, bin_Trt_48K = np.histogram(
        Trt_48K[:, Trt48KIndex], bin_Trt_MC, density=True)
    axes[i].bar(bin_Trt_48K[:-1], hist_Trt_48K * np.diff(bin_Trt_48K)*100, width=np.diff(bin_Trt_48K)[0],
                edgecolor='black', label='48K sims', alpha=0.7)
    axes[i].bar(bin_Trt_MC[:-1], hist_Trt_MC * np.diff(bin_Trt_MC)*100, width=np.diff(bin_Trt_MC)[0],
                edgecolor='black', label='SM, 1M MCs, P='+str(PolyOrder[i].astype(int)), alpha=0.7)
    axes[i].legend()
    #axes[plot_indc[l]].set_xlabel('Bound circulation $\Gamma$ [$m^2/s$]')
    axes[i].set_xlabel('Thrust [kN]', fontsize=20)
    axes[i].grid()
    # print(PCEIndex[j,i])
    l = l+1

#axes[0].set_xlabel('Thrust [kN]', fontsize=20)
#axes[1].set_xlabel('Thrust [kN]', fontsize=20)
axes[0].set_ylabel('Frequency of Occurances [\%]', fontsize=20)



fig3.savefig(timestr+'_MCPCEvsSim_Hist_Trt.pdf', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='pdf', transparent=False)
fig3.savefig(timestr+'_MCPCEvsSim_Hist_Trt.png', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='png', transparent=False)

# %%

f2 = open('20220308_Torque_MC_sims_P2_1M_rlz.pckl', 'rb')
f3 = open('20220308_Torque_MC_sims_P3_1M_rlz.pckl', 'rb')
f4 = open('20220324_Torque_MC_sims_P4_1M_rlz.pckl', 'rb')
f5 = open('NREL5MW_trq_12mps_48000seeds.pckl', 'rb')

trq_MC = np.zeros((4, 1000000))

trq_MC[0, :] = pickle.load(f2)
trq_MC[1, :] = pickle.load(f3)
trq_MC[2, :] = pickle.load(f4)

trq_48K = pickle.load(f5)


f2.close()
f3.close()
f5.close()


trq_MC = trq_MC/1000
trq_48K = trq_48K/1000


mean = np.mean(trq_MC)
standard_deviation = np.std(trq_MC)
distance_from_mean = abs(trq_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
trq_MC_no_ol = trq_MC[not_outlier]


plt.close('all')

plt.rc('font', size=18)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

visib = 0.3
fig3, axes = plt.subplots(1, 3, figsize=(27, 6.35), sharey=True)
timestr = time.strftime('%Y%m%d')

k = 0

plot_indc = ()

for i in range(0, 1):
    for j in range(0, 2):
        k = k+1
        plot_indc = plot_indc+((i, j),)

ITStepSecond = 30
ITStep = ITStepSecond*10
# This is based on the polynomial order and the number of random variables
NoOfRlz = np.array([132, 536, 2002, 6006])
SimLength = np.array([27.3, 10, 10, 10])*10+ITStep
SimLength_xaxis = np.array([27.3, 10, 10, 10])*10
PCEsimLength = np.array([27.3, 6.8, 1.8, 1])
PolyOrder = np.array([2, 3, 4])
elt = np.array([9.942, 40.682, 144.025, 445.721])
l = 0
# PCEIndex[i].astype(int)

trq48KIndex = np.random.randint(60, 5999)

lb = np.min(trq_MC_no_ol)
ub = np.max(trq_MC_no_ol)

bins = np.linspace(lb, ub, 40)


for i in range(0, PolyOrder.size):
    hist_trq_MC, bin_trq_MC = np.histogram(
        trq_MC[i, :], bins=bins, density=True)
    hist_trq_48K, bin_trq_48K = np.histogram(
        trq_48K[:, trq48KIndex], bin_trq_MC, density=True)
    axes[i].bar(bin_trq_48K[:-1], hist_trq_48K * np.diff(bin_trq_48K)*100, width=np.diff(bin_trq_48K)[0],
                edgecolor='black', label='48K sims', alpha=0.7)
    axes[i].bar(bin_trq_MC[:-1], hist_trq_MC * np.diff(bin_trq_MC)*100, width=np.diff(bin_trq_MC)[0],
                edgecolor='black', label='SM, 1M MCs, P='+str(PolyOrder[i].astype(int)), alpha=0.7)
    axes[i].legend()
    #axes[plot_indc[l]].set_xlabel('Bound circulation $\Gamma$ [$m^2/s$]')
    axes[i].set_xlabel('Torque [kNm]', fontsize=20)
    axes[i].grid()
    # print(PCEIndex[j,i])
    l = l+1

#axes[0].set_xlabel('Torque [kNm]', fontsize=20)
axes[0].set_ylabel('Frequency of Occurances [\%]', fontsize=20)



fig3.savefig(timestr+'_MCPCEvsSim_Hist_trq.pdf', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='pdf', transparent=False)
fig3.savefig(timestr+'_MCPCEvsSim_Hist_trq.png', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='png', transparent=False)

# %%

f2 = open('20220308_Thrust_MC_sims_P2_10M_rlz.pckl', 'rb')
f3 = open('20220308_Thrust_MC_sims_P3_10M_rlz.pckl', 'rb')
f4 = open('20220324_Thrust_MC_sims_P4_10M_rlz.pckl', 'rb')

f5 = open('NREL5MW_trt_12mps_48000seeds.pckl', 'rb')

Trt_MC = np.zeros((4, 10000000))

Trt_MC[0, :] = pickle.load(f2)
Trt_MC[1, :] = pickle.load(f3)
Trt_MC[2, :] = pickle.load(f4)


Trt_48K = pickle.load(f5)


f2.close()
f3.close()
f5.close()


Trt_MC = Trt_MC/1000
Trt_48K = Trt_48K/1000


mean = np.mean(Trt_MC)
standard_deviation = np.std(Trt_MC)
distance_from_mean = abs(Trt_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
Trt_MC_no_ol = Trt_MC[not_outlier]


plt.close('all')

plt.rc('font', size=18)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

visib = 0.3
fig3, axes = plt.subplots(3, 1, figsize=(7, 18), sharex=True)
timestr = time.strftime('%Y%m%d')


k = 0

plot_indc = ()

for i in range(0, 1):
    for j in range(0, 2):
        k = k+1
        plot_indc = plot_indc+((i, j),)

ITStepSecond = 30
ITStep = ITStepSecond*10
# This is based on the polynomial order and the number of random variables
NoOfRlz = np.array([132, 536, 2002, 6006])
SimLength = np.array([27.3, 10, 10, 10])*10+ITStep
SimLength_xaxis = np.array([27.3, 10, 10, 10])*10
PCEsimLength = np.array([27.3, 6.8, 1.8, 1])
PolyOrder = np.array([2, 3, 4])
elt = np.array([9.942, 40.682, 144.025, 445.721])
l = 0
# PCEIndex[i].astype(int)

Trt48KIndex = np.random.randint(60, 5999)

lb = np.min(Trt_MC_no_ol)
ub = np.max(Trt_MC_no_ol)

bins = np.linspace(lb, ub, 40)


for i in range(0, PolyOrder.size):
    hist_Trt_MC, bin_Trt_MC = np.histogram(
        Trt_MC[i, :], bins=bins, density=True)
    hist_Trt_48K, bin_Trt_48K = np.histogram(
        Trt_48K[:, Trt48KIndex], bin_Trt_MC, density=True)
    axes[i].bar(bin_Trt_48K[:-1], hist_Trt_48K * np.diff(bin_Trt_48K)*100, width=np.diff(bin_Trt_48K)[0],
                edgecolor='black', label='48K sims', alpha=0.7)
    axes[i].bar(bin_Trt_MC[:-1], hist_Trt_MC * np.diff(bin_Trt_MC)*100, width=np.diff(bin_Trt_MC)[0],
                edgecolor='black', label='SM, 10M MCs, P='+str(PolyOrder[i].astype(int)), alpha=0.7)
    axes[i].legend()
    #axes[plot_indc[l]].set_xlabel('Bound circulation $\Gamma$ [$m^2/s$]')
    axes[i].set_ylabel('Frequency of Occurances [\%]', fontsize=20)
    axes[i].grid()
    # print(PCEIndex[j,i])
    l = l+1

#axes[0].set_xlabel('Thrust [kN]', fontsize=20)
axes[2].set_xlabel('Thrust [kN]', fontsize=20)


fig3.savefig(timestr+'_10M_MCPCEvsSim_Hist_Trt.pdf', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='pdf', transparent=False)
fig3.savefig(timestr+'_10M_MCPCEvsSim_Hist_Trt.png', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='png', transparent=False)

# %%

f2 = open('20220308_Torque_MC_sims_P2_10M_rlz.pckl', 'rb')
f3 = open('20220308_Torque_MC_sims_P3_10M_rlz.pckl', 'rb')
f4 = open('20220324_Torque_MC_sims_P4_10M_rlz.pckl', 'rb')
f5 = open('NREL5MW_trq_12mps_48000seeds.pckl', 'rb')

trq_MC = np.zeros((4, 10000000))

trq_MC[0, :] = pickle.load(f2)
trq_MC[1, :] = pickle.load(f3)
trq_MC[2, :] = pickle.load(f4)

trq_48K = pickle.load(f5)


f2.close()
f3.close()
f5.close()


trq_MC = trq_MC/1000
trq_48K = trq_48K/1000


mean = np.mean(trq_MC)
standard_deviation = np.std(trq_MC)
distance_from_mean = abs(trq_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
trq_MC_no_ol = trq_MC[not_outlier]


plt.close('all')

plt.rc('font', size=18)
plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18)

visib = 0.3
fig3, axes = plt.subplots(3, 1, figsize=(7, 18), sharex=True)
timestr = time.strftime('%Y%m%d')


k = 0

plot_indc = ()

for i in range(0, 1):
    for j in range(0, 2):
        k = k+1
        plot_indc = plot_indc+((i, j),)

ITStepSecond = 30
ITStep = ITStepSecond*10
# This is based on the polynomial order and the number of random variables
NoOfRlz = np.array([132, 536, 2002, 6006])
SimLength = np.array([27.3, 10, 10, 10])*10+ITStep
SimLength_xaxis = np.array([27.3, 10, 10, 10])*10
PCEsimLength = np.array([27.3, 6.8, 1.8, 1])
PolyOrder = np.array([2, 3, 4])
elt = np.array([9.942, 40.682, 144.025, 445.721])
l = 0
# PCEIndex[i].astype(int)

trq48KIndex = np.random.randint(60, 5999)

lb = np.min(trq_MC_no_ol)
ub = np.max(trq_MC_no_ol)

bins = np.linspace(lb, ub, 40)


for i in range(0, PolyOrder.size):
    hist_trq_MC, bin_trq_MC = np.histogram(
        trq_MC[i, :], bins=bins, density=True)
    hist_trq_48K, bin_trq_48K = np.histogram(
        trq_48K[:, trq48KIndex], bin_trq_MC, density=True)
    axes[i].bar(bin_trq_48K[:-1], hist_trq_48K * np.diff(bin_trq_48K)*100, width=np.diff(bin_trq_48K)[0],
                edgecolor='black', label='48K sims', alpha=0.7)
    axes[i].bar(bin_trq_MC[:-1], hist_trq_MC * np.diff(bin_trq_MC)*100, width=np.diff(bin_trq_MC)[0],
                edgecolor='black', label='SM, 10M MCs, P='+str(PolyOrder[i].astype(int)), alpha=0.7)
    axes[i].legend()
    #axes[plot_indc[l]].set_xlabel('Bound circulation $\Gamma$ [$m^2/s$]')
    axes[i].set_ylabel('Frequency of Occurances [\%]', fontsize=20)
    axes[i].grid()
    # print(PCEIndex[j,i])
    l = l+1

#axes[0].set_xlabel('Torque [kNm]', fontsize=20)
axes[2].set_xlabel('Torque [kNm]', fontsize=20)


fig3.savefig(timestr+'_10M_MCPCEvsSim_Hist_trq.pdf', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='pdf', transparent=False)
fig3.savefig(timestr+'_10M_MCPCEvsSim_Hist_trq.png', dpi=300, facecolor='w',
             edgecolor='w', orientation='portrait', format='png', transparent=False)

# %%

f2 = open('20220308_Thrust_MC_sims_P2_1M_rlz.pckl', 'rb')
f3 = open('20220308_Thrust_MC_sims_P3_1M_rlz.pckl', 'rb')
f4 = open('20220324_Thrust_MC_sims_P4_1M_rlz.pckl', 'rb')

f5 = open('NREL5MW_trt_12mps_48000seeds.pckl', 'rb')

Trt_MC = np.zeros((3, 1000000))

Trt_MC[0, :] = pickle.load(f2)
Trt_MC[1, :] = pickle.load(f3)
Trt_MC[2, :] = pickle.load(f4)




Trt_48K = pickle.load(f5)

f2.close()
f3.close()
f4.close()
f5.close()


Ntimes = Trt_48K.shape[1]

mean = np.mean(Trt_MC)
standard_deviation = np.std(Trt_MC)
distance_from_mean = abs(Trt_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
Trt_MC_no_ol = Trt_MC[not_outlier]



lb = np.min(Trt_MC_no_ol)
ub = np.max(Trt_MC_no_ol)

bins_o = np.linspace(lb,ub,30)

Trt_hist = np.zeros((bins_o.shape[0]-1,Ntimes))
Trt_MC_hist = np.zeros((bins_o.shape[0]-1,Trt_MC.shape[0]))
HL_dist= np.zeros((Trt_MC.shape[0],Ntimes))
Trt_hist = np.zeros(bins_o.shape[0]-1)
Trt_MC_hist = np.zeros((bins_o.shape[0]-1,Trt_MC.shape[0]))
#%%
for i in range(0,Trt_MC.shape[0]):
    #lb = np.min()
    #ub = np.max(Trt_MC[i,:])
    bins_i = np.linspace(lb,ub,30)
    Trt_MC_hist[:,i],_ = np.histogram(Trt_MC[i,:],bins=bins_i,density='True')
    Q = np.sqrt(Trt_MC_hist*np.diff(bins_o)[0])
    #Q = Trt_MC_hist
    for j in range(0,Ntimes):    
        Trt_hist,blabla = np.histogram(Trt_48K[:,j],bins=bins_i,density='True')
        P = np.sqrt(Trt_hist*np.diff(bins_o)[0])
        #P = Trt_hist
        HH = (P - Q[:,i])**2
        HL_dist[i,j] = (1/np.sqrt(2))*(np.sqrt(np.sum(HH)))
        #HL_dist[i,j] = st.entropy(P, Q[:,i])
        print(i,j)

timestr = time.strftime('%Y%m%d')  

f=open(timestr+'_HellingerDist_Trt_vs_SQGMC.pckl','wb')
pickle.dump(HL_dist,f)
f.close()    

#%%

f=open('20220425_HellingerDist_Trt_vs_SQGMC.pckl','rb')
HL_dist = pickle.load(f)
f.close()    

timestr = time.strftime('%Y%m%d')  
plt.close('all')
t_fft = np.linspace(0,600,Ntimes)

plt.close('all')    

visib=0.3
fig, axes = plt.subplots(1,1,figsize=(11.69,8.27),sharey=True)
PolyOrder = np.array([2,3,4])


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

fig.savefig(timestr+'_HL_Dis_time_Trt_vs_SQGMCs.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr+'_HL_Dis_time_Trt_vs_SQGMCs.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)

#%%

f2 = open('20220308_Torque_MC_sims_P2_1M_rlz.pckl', 'rb')
f3 = open('20220308_Torque_MC_sims_P3_1M_rlz.pckl', 'rb')
f4 = open('20220324_Torque_MC_sims_P4_1M_rlz.pckl', 'rb')

f5 = open('NREL5MW_trq_12mps_48000seeds.pckl', 'rb')

Trq_MC = np.zeros((3, 1000000))

Trq_MC[0, :] = pickle.load(f2)
Trq_MC[1, :] = pickle.load(f3)
Trq_MC[2, :] = pickle.load(f4)




Trq_48K = pickle.load(f5)

f2.close()
f3.close()
f4.close()
f5.close()


Ntimes = Trq_48K.shape[1]

mean = np.mean(Trq_MC)
standard_deviation = np.std(Trq_MC)
distance_from_mean = abs(Trq_MC - mean)
max_deviations = 4.5
not_outlier = distance_from_mean < max_deviations * standard_deviation
Trq_MC_no_ol = Trq_MC[not_outlier]



lb = np.min(Trq_MC_no_ol)
ub = np.max(Trq_MC_no_ol)

bins_o = np.linspace(lb,ub,30)

Trq_hist = np.zeros((bins_o.shape[0]-1,Ntimes))
Trq_MC_hist = np.zeros((bins_o.shape[0]-1,Trq_MC.shape[0]))
HL_dist= np.zeros((Trq_MC.shape[0],Ntimes))
Trq_hist = np.zeros(bins_o.shape[0]-1)
Trq_MC_hist = np.zeros((bins_o.shape[0]-1,Trq_MC.shape[0]))

for i in range(0,Trq_MC.shape[0]):
    #lb = np.min()
    #ub = np.max(Trq_MC[i,:])
    bins_i = np.linspace(lb,ub,30)
    Trq_MC_hist[:,i],_ = np.histogram(Trq_MC[i,:],bins=bins_i,density='True')
    Q = np.sqrt(Trq_MC_hist*np.diff(bins_o)[0])
    #Q = Trq_MC_hist
    for j in range(0,Ntimes):    
        Trq_hist,blabla = np.histogram(Trq_48K[:,j],bins=bins_i,density='True')
        P = np.sqrt(Trq_hist*np.diff(bins_o)[0])
        #P = Trq_hist
        HH = (P - Q[:,i])**2
        HL_dist[i,j] = (1/np.sqrt(2))*(np.sqrt(np.sum(HH)))
        #HL_dist[i,j] = st.entropy(P, Q[:,i])
        print(i,j)

timestr = time.strftime('%Y%m%d')  

f=open(timestr+'_HellingerDist_Trq_vs_SQGMC.pckl','wb')
pickle.dump(HL_dist,f)
f.close()    

timestr = time.strftime('%Y%m%d')  
plt.close('all')
t_fft = np.linspace(0,600,Ntimes)

plt.close('all')    

visib=0.3
fig, axes = plt.subplots(1,1,figsize=(11.69,8.27),sharey=True)
PolyOrder = np.array([2,3,4])


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

fig.savefig(timestr+'_HL_Dis_time_Trq_vs_SQGMCs.pdf',dpi=300,facecolor='w', 
             edgecolor='w',orientation='portrait', format='pdf',transparent=False)
fig.savefig(timestr+'_HL_Dis_time_Trq_vs_SQGMCs.png',dpi=300,facecolor='w', 
                 edgecolor='w',orientation='portrait', format='png',transparent=False)



