#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:19:58 2022

@author: noahbungart
"""

#Import packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
# %matplotlib inline
import plotly.express as px
import plotly.io as pio

#Import COSMIC and other things
import cosmic
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve
from cosmic.sample.sampler import independent

##############################################################################

#Defining function for rv_variability and filtering
def rv_variable(m1, m2, a, period, ecc, sin_i):
    """
    Function to calculate readial velocity variability
    
    m1: Mass 1
    m2: Mass 2
    period: Period
    ecc: Eccentricity
    a: amplitude
    """
    var = (2*np.pi*a*m2*sin_i)/(period*(m1+m2)*(1-ecc**2)**(1/2))
    return var

def filter_func(m1 , m2, period , ecc, a):
    """
    Function that filters out unwanted values for each parameter 
    Such as 0,-1,inf,nan
    """
    p = period[period != 0]
    p = p[p != -1]
    p = p[p != np.inf]
    
    semi = a[a != 0]
    semi = semi[semi != -1]
    
    e = ecc[ecc != -1]
    
    #Start the filtering!
    #Can check the del_arr as it goes by uncommenting the print lines
    
    #period indecies
    x = period.index[period == 0]
    y = period.index[period == -1]
    z = period.index[period == np.inf]

    #Update del_arr
    del_arr = x
    del_arr = del_arr.append(y)
    del_arr = del_arr.append(z)
    #print(del_arr)

    #Semi major indecies
    x_2 = a.index[a == 0]
    y_2 = a.index[a == -1]

    #Update del_arr
    del_arr = del_arr.append(x_2)
    del_arr = del_arr.append(y_2)
    #print(del_arr)

    #Ecc indecies
    x_3 = ecc.index[ecc == -1]

    #Update del_arr
    del_arr = del_arr.append(x_3)
    #print(del_arr)
    
    #Create final array and remove duplicates
    delete_arr = np.unique(del_arr)

    ecc_f = []
    for i in range(len(ecc)):
        if ecc.index[i] not in delete_arr:
            ecc_f.append(ecc[ecc.index[i]])
    period_f = []
    for i in range(len(period)):
        if period.index[i] not in delete_arr:
            period_f.append(period[period.index[i]])
    a_f = []
    for i in range(len(a)):
        if a.index[i] not in delete_arr:
            a_f.append(a[a.index[i]])
    
    #Update the masses
    m1 = []
    for i in range(len(mass1)):
        if mass1.index[i] not in delete_arr:
            m1.append(mass1[mass1.index[i]])
    m2 = []
    for i in range(len(mass2)):
        if mass2.index[i] not in delete_arr:
            m2.append(mass2[mass2.index[i]])
            
    return m1, m2, period_f, ecc_f, a_f

##############################################################################
#Set the runID, what change is being made, and size
mainID = 'Ratio Plots'
param_change = 'porb_model'
size = 100000
runID = 'log_uniform'


#Check paths and make new directories if necessary
path1 = mainID + '/'
isExist = os.path.exists(path1)

if not isExist:
  
  # Create a new directory because it does not exist 
  os.makedirs(path1)
  print("The new directory is created!")

path2 = path1 + param_change + '/'
isExist = os.path.exists(path2)

if not isExist:
  
  # Create a new directory because it does not exist 
  os.makedirs(path2)
  print("The new directory is created!")

path3 = path2 + str(size) + '/'
isExist = os.path.exists(path3)

if not isExist:
  
  # Create a new directory because it does not exist 
  os.makedirs(path3)
  print("The new directory is created!")

path4 = path3 + runID + '/'
isExist = os.path.exists(path4)

if not isExist:
  
  # Create a new directory because it does not exist 
  os.makedirs(path4)
  print("The new directory is created!")

##############################################################################
#Set Initial binary parameters

#Setting what evolution types are allowed
final_kstar1 = [10,11,12]
final_kstar2 = [10,11,12]

#Set the initial binary population parameters
InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = \
     InitialBinaryTable.sampler('independent', final_kstar1, final_kstar2, binfrac_model=1.,
                                primary_model='kroupa01', ecc_model='sana12', porb_model='log_uniform',
                                qmin=-1, SF_start=13700.0, SF_duration=0.0, met=0.02, size=size)

#Can print initial binaries here to check
#print(InitialBinaries)

print('####################### Initial Binaries Set ##########################')

#Set the BSEDict
BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'dtp' : 13700.0}

#Evolve the system
bpp, bcm, initC, kick_info  = Evolve.evolve(initialbinarytable=InitialBinaries, BSEDict=BSEDict)

print(bcm.iloc[:10])

print('###################### System Evolved ################################')

#Get all parameters and create sini artificial data
mass1 = bcm.mass_1[bcm.tphys == 13700.0]
#print(mass1)
mass2 = bcm.mass_2[bcm.tphys == 13700.0]
#print(mass2)
period = bcm.porb[bcm.tphys == 13700.0]
#print(period)
ecc = bcm.ecc[bcm.tphys == 13700.0]
#print(ecc)
a = bcm.sep[bcm.tphys == 13700.0]
#print(a)

#Checking to make sure all initial dataframes are the 
#same length
#print(len(mass1))
#print(len(mass2))
#print(len(period))
#print(len(ecc))
#print(len(a))

m1, m2, period_f, ecc_f, a_f = filter_func(mass1, mass2, period, ecc, a)

sini = np.random.uniform(0, 1, len(m1))

print('######################### Filtering Complete #######################')

#Make the dataframes arrays in order to use rv funciton
m1 = np.array(m1)
m2 = np.array(m2)
ecc_f = np.array(ecc_f)
period_f = np.array(period_f)
a_f = np.array(a_f)
sini = np.array(sini)

rv = rv_variable(m1, m2, a_f, period_f, ecc_f, sini)

rv = rv[~np.isnan(rv)]

#Can check rv's below if wanted
#print(rv)
#print(len(rv))

print('########################### Start Plotting ##########################')

q = m1 / m2

m1 = m1[0:len(rv)]
a_f = a_f[0:len(rv)]
period_f = period_f[0:len(rv)]
q = q[0:len(rv)]

fig1 = px.scatter(x=m1, y=q, color=rv, range_color=[0,20],
                  labels={'x' :'Mass 1', 'y' :'Mass Ratio'},
                  title='Mass Ratio vs M1')
fig1.write_image(path4 + 'QvsM1_' + runID +'.png')

fig2 = px.scatter(x=a_f, y=q, color=rv, range_color=[0,10], log_x=True,
                  range_x = [10000,500000],
                  range_y = [0, 3],
                  labels={'x' :'Semi Major', 'y' :'Mass Ratio'},
                  title='Mass Ratio vs Semi-major')
fig2.write_image(path4 + 'QvsSemi_' + runID +'.png')

fig3 = px.scatter(x=period_f, y=q, color=rv, range_color=[0,10],log_x=True,
                  range_x = [1000000,30000000],
                  range_y = [0, 3],
                  labels={'x':'Period', 'y':'Mass Ratio'},
                  title='Mass Ratio vs Period')
fig3.write_image(path4 + 'QvsPeriod_' + runID +'.png')

print('######################### Get largest Mass Ratio ####################')
mass1 = bcm.mass_1[bcm.tphys == 13700.0]
mass2 = bcm.mass_2[bcm.tphys == 13700.0]
q = mass1 / mass2

q[np.isinf(q)] = 0

max_ind = np.argmax(q)
print("Largest mass ratio:", q[max_ind])
print(bcm.merger_type[bcm.tphys == 13700.0][max_ind])
print(bcm.SN_1[bcm.tphys == 13700.0][max_ind])
print(bcm.bin_state[bcm.tphys == 13700.0][max_ind])