#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 20 15:41:22 2021
This file is to run BEM simulations in batches. Hopefully it works.
In this file I run simulations to fit PCE to the output.
This is for 6006 simulations.
@author: rad
"""

import chaospy as cp
import numpy as np
import pickle
from datetime import date
import time
import bemused
from makewindVeers_main import *
import os

#%% Unsteady wind, one set to make the matrices


Nb = 1
dt = 0.1;
rho = 1.225
plots =0
x=0
y=0

v0=8
Nf=10
height = 100

xi = np.random.random((Nb,Nf))
P = np.vstack((x,y))

v0_unsteady,t_fft,Uk,_,_,_,_,Fc,X =  makeWindVeers_630s (v0,P,Nf,dt,plots,height,xi)
v0_unsteady = v0_unsteady.T


#%% NREL 5MW details and model

blade_definition = bemused.Blade.from_yaml('NREL5MW_simpblade.yaml')
aerofoil_database_filename = 'oc3_aerofoils.npz'
aerofoil_database = bemused.AerofoilDatabase(aerofoil_database_filename)
root_length = 1.5

# Let's have a normal three-bladed rotor.
number_of_blades = 3

model = bemused.BEMModel(
  blade_definition, root_length, number_of_blades, aerofoil_database)


f = open('NREL_5MW_Data.pckl','rb')
NREL5MW_Data = pickle.load(f)
f.close()

wind_speeds = NREL5MW_Data[:,0]
pitch_angle = NREL5MW_Data[:,7]*(np.pi/180)
rotor_speed = NREL5MW_Data[:,6]*(np.pi/30)



thrust_wspd = np.zeros((len(wind_speeds),len(t_fft)))
torque_wspd = np.zeros((len(wind_speeds),len(t_fft)))
un_wspd = np.zeros((len(wind_speeds),len(t_fft)))

#%%

NoOfRlz = np.array([48000]) #This is based on the polynomial order and the number of random variables


for Seeds in NoOfRlz:
    xi = cp.generate_samples(Seeds,domain=cp.Iid(cp.Uniform(0,1),10),rule = 'S')
    newpath = (os.getcwd() + '/NREL5MW_BEM_UnSteady_NoOfSims_' +str(Seeds))
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    file_path_xi = newpath+'/xi_'+str(Seeds)+'seeds.pckl'
    f = open(file_path_xi ,'wb')
    pickle.dump(xi,f)
    f.close()
    for i in range(0,len(wind_speeds)):        
        BEM_FrozenWake = bemused.models.FrozenWakeAerodynamics(model, wind_speeds[i],rotor_speed[i], pitch_angle[i])
        thrust_wspd = np.zeros((Seeds,len(t_fft)))
        torque_wspd = np.zeros((Seeds,len(t_fft)))
        for j in range(0,Seeds):
            forces = np.zeros((v0_unsteady.shape[0], len(blade_definition.x), 2))
            v0_unsteady,t_fft,Uk,_,_,_,_,Fc,X =  makeWindVeers_630s (wind_speeds[i],P,Nf,dt,0,height,xi[:,j].reshape(1,-1))
            v0_unsteady = v0_unsteady.flatten()
            for k in range(0,v0_unsteady.shape[0]):
                forces[k] =  BEM_FrozenWake.forces(v0_unsteady[k], rotor_speed[i], pitch_angle[i],rho)
            thrust_wspd[j,:] = number_of_blades * np.trapz(forces[:, :, 0], model.radii, axis=1)
            torque_wspd[j,:] = (-number_of_blades *np.trapz(model.radii * forces[:, :, 1], model.radii, axis=1))
            print(wind_speeds[i],j)
        file_path_torque = newpath+'/NREL5MW_trq_'+str(wind_speeds[i].astype(int))+'mps_'+str(Seeds)+'seeds.pckl'
        file_path_thrust = newpath+'/NREL5MW_trt_'+str(wind_speeds[i].astype(int))+'mps_'+str(Seeds)+'seeds.pckl'
        f = open(file_path_torque,'wb')
        pickle.dump(torque_wspd,f)
        f.close()
        f = open(file_path_thrust,'wb')
        pickle.dump(thrust_wspd,f)
        f.close()