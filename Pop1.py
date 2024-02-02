# -*- coding: utf-8 -*-
"""
Population Editing

This is where I modify population parameters and see how it affects the
distribution of WD-WD binaries
"""
#Import necessary packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import cosmic
from cosmic.sample.initialbinarytable import InitialBinaryTable
from cosmic.evolve import Evolve
from cosmic.sample.sampler import independent

#Set evolution types to only be WD-WD or merge
final_kstar1 = [10,11,12]
final_kstar2 = [10,11,12]

#Set the binary parameters
InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', 
                                                                                                 final_kstar1, 
                                                                                                 final_kstar2, 
                                                                                                 binfrac_model=0.5, 
                                                                                                 primary_model='kroupa01', 
                                                                                                 ecc_model='sana12', 
                                                                                                 porb_model='sana12', 
                                                                                                 qmin=-1, 
                                                                                                 SF_start=13700.0, 
                                                                                                 SF_duration=0.0, 
                                                                                                 met=0.02, 
                                                                                                 size=100000)
print(InitialBinaries)

BSEDict = {'xi': 1.0, 'bhflag': 1, 'neta': 0.5, 'windflag': 3, 'wdflag': 1, 'alpha1': 1.0, 'pts1': 0.001, 'pts3': 0.02, 'pts2': 0.01, 'epsnov': 0.001, 'hewind': 0.5, 'ck': 1000, 'bwind': 0.0, 'lambdaf': 0.0, 'mxns': 3.0, 'beta': -1.0, 'tflag': 1, 'acc2': 1.5, 'grflag' : 1, 'remnantflag': 4, 'ceflag': 0, 'eddfac': 1.0, 'ifflag': 0, 'bconst': 3000, 'sigma': 265.0, 'gamma': -2.0, 'pisn': 45.0, 'natal_kick_array' : [[-100.0,-100.0,-100.0,-100.0,0.0], [-100.0,-100.0,-100.0,-100.0,0.0]], 'bhsigmafrac' : 1.0, 'polar_kick_angle' : 90, 'qcrit_array' : [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'cekickflag' : 2, 'cehestarflag' : 0, 'cemergeflag' : 0, 'ecsn' : 2.25, 'ecsn_mlow' : 1.6, 'aic' : 1, 'ussn' : 0, 'sigmadiv' :-20.0, 'qcflag' : 1, 'eddlimflag' : 0, 'fprimc_array' : [2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0,2.0/21.0], 'bhspinflag' : 0, 'bhspinmag' : 0.0, 'rejuv_fac' : 1.0, 'rejuvflag' : 0, 'htpmb' : 1, 'ST_cr' : 1, 'ST_tide' : 1, 'bdecayfac' : 1, 'rembar_massloss' : 0.5, 'kickflag' : 0, 'zsun' : 0.014, 'bhms_coll_flag' : 0, 'don_lim' : -1, 'acc_lim' : -1, 'dtp' : 13700.0}

bpp, bcm, initC, kick_info = Evolve.evolve(initialbinarytable = InitialBinaries, BSEDict=BSEDict)

print(bcm.iloc[:10])

#print(bpp)
#print(bcm)

#Making our arrays for necessary data
mass1=bcm.mass_1
m1=np.array(mass1)
mass2=bcm.mass_2
m2=np.array(mass2)

tphys=bcm.tphys
Tphys=np.array(tphys)

final_mass1 = [m1[i] for i in range(len(Tphys)) if Tphys[i] == 13700.0]
f_m1 = pd.DataFrame(final_mass1)
final_mass2 = [m2[i] for i in range(len(Tphys)) if Tphys[i] == 13700.0]
f_m2 = pd.DataFrame(final_mass2)
#print(final_mass1)
#print(final_mass2)

#Start plotting the histograms
f_m1.plot(kind='hist', edgecolor='black', alpha=1, legend=None)
plt.title('Final Masses of Star 1')
plt.xlabel('Range of Masses')
plt.show()

f_m2.plot(kind='hist', edgecolor='black', alpha=1, legend=None)
plt.title('Final Masses of Star 2')
plt.xlabel('Range of Masses')
plt.savefig("Figs/Star 2 Masses")
plt.show()

#Removing binary systems with period of infinity
fperiod = bcm.porb
fp = np.array(fperiod)
f_period = [fp[i] for i in range(len(Tphys)) if Tphys[i] == 13700.0]
f = [f_period[i] for i in range(len(f_period)) if f_period[i]!= np.Inf]
fp1 = pd.DataFrame(f)

fp1.plot(kind='hist', edgecolor='black', alpha=1, legend=None)
plt.title('Final Period of Stars')
plt.xlabel('Periods')

print("Complete :)")