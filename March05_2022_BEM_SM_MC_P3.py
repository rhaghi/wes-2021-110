#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 20 13:27:58 2021

@author: rad
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 17:19:41 2021
This is to fit SM on the sectional statsitics for the WES paper.
@author: rad
"""

#
import chaospy as cp
import numpy as np
#import matplotlib.pyplot as plt
import pickle
#from makewindVeers_main import *
#import seaborn as sns
#import scipy.stats as st
#import statistics as sts
#import scipy.fftpack
#import scipy.special
#import scipy.io as sio
from datetime import date
import time
#from mpl_toolkits.mplot3d import Axes3D
#from mpl_toolkits.mplot3d.art3d import Poly3DCollection # New import
#import pandas as pd
#import rainflow
#import fatpack
import os
#import sklearn.linear_model as sklm
import sys



#plt.rc('text', usetex=True)
#plt.rc('font', family='serif')
#plt.rc('font', size = 14)
#plt.rc('xtick', labelsize = 14)
#%%

#%%

#timestr = '20200622'
ITStepSecond = 30
ITStep = ITStepSecond*10
#NoOfRlz = np.array([22,132,536])#,2002,6006]) #This is based on the polynomial order and the number of random variables
#NoOfRlz = np.array([2002]) #To test the code
#SimLength = np.array([163.4,27.3,10,10,10])
SimLength = np.array([10])#,10,10,10])
PCEIndex = np.round(SimLength/2)
#SimLength = np.array([0.1]) #To test the code
PCEIndex= PCEIndex+ITStepSecond 
PCEIndex = PCEIndex *10
polyOrder = np.array([3])#,3,4,5])
#polyOrder = np.array([1])#To test the code
MC_sims_no = np.array([10,100,1000,10000,48000,1e5,5e5,1e6,1e7,1e8])
dist =InitialDist = cp.Iid(cp.Uniform(0,1),10)
cpu_time_trt= np.zeros((MC_sims_no.size,polyOrder.size))
cpu_time_trq= np.zeros((MC_sims_no.size,polyOrder.size))
poly_load_thrust = cp.load('BEM_Thrust_SM_PolyOrder_3.npz',allow_pickle=True)
poly_load_torque = cp.load('BEM_Torque_SM_PolyOrder_3.npz',allow_pickle=True)

PCE_trt = poly_load_thrust['arr_0'].item()
PCE_trq = poly_load_torque['arr_0'].item()

#PCEIndex = np.round(SimLength/2)

tstr_file = time.strftime('%Y%m%d')

for i in range(0,MC_sims_no.size):
    for j in range(0,polyOrder.size):
        #print(MC_sims_no[i],polyOrder[j])
        t0 = time.time()
        Thrust_MC = cp.call(PCE_trt[str(j),str(PCEIndex[j].astype(int))],cp.generate_samples(MC_sims_no[i].astype(int),InitialDist))
        cpu_time_trt[i,j]=time.time()-t0
        Torque_MC = cp.call(PCE_trq[str(j),str(PCEIndex[j].astype(int))],cp.generate_samples(MC_sims_no[i].astype(int),InitialDist))
        cpu_time_trq[i,j]=time.time()-t0
        if i>6:
            f = open(tstr_file+'_Thrust_MC_sims_P'+str(polyOrder[j].astype(int))+'_'+str((MC_sims_no[i]/1000000).astype(int))+'M_rlz.pckl','wb')
            pickle.dump(Thrust_MC,f)
            f.close()
            f = open(tstr_file+'_Torque_MC_sims_P'+str(polyOrder[j].astype(int))+'_'+str((MC_sims_no[i]/1000000).astype(int))+'M_rlz.pckl','wb')
            pickle.dump(Torque_MC,f)
            f.close()
        sys.stdout.write(str(MC_sims_no[i])+'  '+str(polyOrder[j])+'  ' +'\n')
        sys.stdout.flush()

f=open(tstr_file+'_CPU_Time_P3_MC_trt.pckl','wb')
pickle.dump(cpu_time_trt,f)
f.close()

f=open(tstr_file+'_CPU_Time_P3_MC_trq.pckl','wb')
pickle.dump(cpu_time_trq,f)
f.close()

