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
wind_speed = 12

#timestr = '20200622'
ITStepSecond = 30
ITStep = ITStepSecond*10
NoOfRlz = np.array([2002]) #This is based on the polynomial order and the number of random variables
#NoOfRlz = np.array([2002]) #To test the code
SimLength = np.array([10])
#SimLength = np.array([0.1]) #To test the code
SimLength= SimLength+ITStepSecond 
SimLength = SimLength *10
polyOrder = np.array([4])
#polyOrder = np.array([4])#To test the code


TrtMeanPCE_ShortSim = np.empty((SimLength.shape[0],int(SimLength.max())))
TrtMeanPCE_ShortSim[:] = None

TrtStdPCE_ShortSim = np.empty((SimLength.shape[0],int(SimLength.max())))
TrtStdPCE_ShortSim[:] = None

TrqMeanPCE_ShortSim = np.empty((SimLength.shape[0],int(SimLength.max())))
TrqMeanPCE_ShortSim[:] = None

TrqStdPCE_ShortSim = np.empty((SimLength.shape[0],int(SimLength.max())))
TrqStdPCE_ShortSim[:] = None


print('setting the poly order')




ApproxAll_trt = {}
ApproxAll_trq = {}

print('about to start the loop')
InitialDist = cp.Iid(cp.Uniform(0,1),10)
#orthPoly = cp.generate_expansion(polyOrder[ii], InitialDist)

cpu_time = np.zeros((len(NoOfRlz),2))

for i in range(0, NoOfRlz.size):
    orthPoly = cp.generate_expansion(polyOrder[i], InitialDist)
    file_path_trt = 'NREL5MW_trt_'+str(wind_speed)+'mps_'+str(NoOfRlz[i].astype(int))+'seeds.pckl'
    file_path_trq = 'NREL5MW_trq_'+str(wind_speed)+'mps_'+str(NoOfRlz[i].astype(int))+'seeds.pckl'
    file_path_xi =  'xi_' + str(NoOfRlz[i].astype(int))+'seeds.pckl'
    f = open(file_path_trt, 'rb')
    Thrust = pickle.load(f)
    f.close()
    f = open(file_path_trq, 'rb')
    Torque = pickle.load(f)
    f.close()
    f = open(file_path_xi, 'rb')
    xi = pickle.load(f)
    f.close()
    cpu_time_loop = np.zeros((SimLength[i].astype(int)-ITStep,2))
    k=0;
    for j in range(ITStep, SimLength[i].astype(int)):
        dataPoints = np.zeros(NoOfRlz[i].astype(int))
        samples_trt = np.zeros(NoOfRlz[i].astype(int))
        samples_trt = Thrust[:, j]
        samples_trq = np.zeros(NoOfRlz[i].astype(int))
        samples_trq = Torque[:, j]
        dataPoints = xi
        t0 = time.time()
        approx_trt = cp.fit_regression(orthPoly, dataPoints, samples_trt)
        cpu_time_loop[k,0] = time.time()-t0
        t0 = time.time()
        approx_trq = cp.fit_regression(orthPoly, dataPoints, samples_trq)
        cpu_time_loop[k,1] = time.time()-t0
        #TrtMeanPCE_ShortSim[i,j] = cp.E(approx_trt,InitialDist)
        #TrtStdPCE_ShortSim[i,j] = cp.Std(approx_trt,InitialDist)
        #TrqMeanPCE_ShortSim[i,j] = cp.E(approx_trq,InitialDist)
        #TrqStdPCE_ShortSim[i,j] = cp.Std(approx_trq,InitialDist)
        ApproxAll_trt[str(i), str(j)] = approx_trt
        ApproxAll_trq[str(i), str(j)] = approx_trq
        sys.stdout.write(str(i)+'  '+str(j)+'  ' +
                         'PolyOrder: '+str(polyOrder[i])+'\n')
        sys.stdout.flush()
        k=k+1
    cp.savez('BEM_Thrust_SM_PolyOrder_'+str(polyOrder[i]), ApproxAll_trt)
    cp.savez('BEM_Torque_SM_PolyOrder_'+str(polyOrder[i]), ApproxAll_trq)
    cpu_time[i,:] = np.mean(cpu_time_loop,0)

f = open('cpu_time_SM_P4.pckl','wb')
pickle.dump(cpu_time,f)
f.close()