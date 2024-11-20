# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 20:12:48 2024

@author: danab
"""

import Microbial_CORPSE_solvers
from numpy import array,arange
import copy

# This makes an array of all the model time steps (in units of years). In this case, it starts at zero, ends at 120 days, and has a time step of one day
t=arange(0,70/365,1/365)

# Parameters controlling the model (sandy soil, no burn)
params={
    #'vmaxref':{'Fast':8.0,'Slow':0.2,'Necro':4.0, 'Py':0.1}, #  Relative maximum enzymatic decomp rates for each C type (year-1)
    'vmaxref':{'MBC_1':{'Fast':3.0,'Slow':0.15,'Necro':7.0, 'Py':0.1},
               'MBC_2':{'Fast':0.0,'Slow':0.0,'Necro':0.0, 'Py':0.0},
               'MBC_3':{'Fast':0.0,'Slow':0.0,'Necro':0.0, 'Py':0.0}, 
               'MBC_4':{'Fast':0.0,'Slow':0.0,'Necro':0.0, 'Py':0.0}},
    'Ea':{'Fast':5e3,'Slow':30e3,'Necro':5e3, 'Py':35e3},      # Activation energy (controls T dependence)
    'kC':{'MBC_1':{'Fast':0.02,'Slow':0.01,'Necro':0.02, 'Py':0.1},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
          'MBC_2':{'Fast':0.02,'Slow':0.1,'Necro':0.02, 'Py':0.1},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
          'MBC_3':{'Fast':0.1,'Slow':0.01,'Necro':0.1, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
          'MBC_4':{'Fast':0.1,'Slow':0.01,'Necro':0.1, 'Py':0.01}},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    # Smaller kC results in more CO2 produced for a given amount of C available to decompose?
    # Copiotrophs should have larger kC than oligotrophs. Oligotrophs 
    #   have higher substrate affinity.
    'gas_diffusion_exp':0.6,  # Determines suppression of decomposition at high soil moisture
    'substrate_diffusion_exp':1.5,   # Controls suppression of decomp at low soil moisture
    'minMicrobeC':{'MBC_1':1e-4,'MBC_2':1e-5,'MBC_3':1e-3,'MBC_4':1e-3}, # Minimum microbial biomass (fraction of total C). Prevents microbes from going extinct and allows some slow decomposition under adverse conditions
    'Tmic':{'MBC_1':0.5,'MBC_2':0.06,'MBC_3':0.25,'MBC_4':0.25}, # Microbial biomass mean lifetime (years)
            # Give MBC 1 a mean lifetime of 6 hours (0.25) and MBC_2 a mean lifetime of 60 minutes.
    'et':{'MBC_1':0.8,'MBC_2':0.8,'MBC_3':0.6,'MBC_4':0.6},                  # Fraction of microbial biomass turnover (death) that goes to necromass instead of being immediately mineralized to CO2
    # Wouldn't we expect all the MBC to become necromass and then decomposition of necromass --> CO2??
    'eup':{'MBC_1':{'Fast':0.5,'Slow':0.4,'Necro':0.6, 'Py':0.1},     # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
           'MBC_2':{'Fast':0.4,'Slow':0.4,'Necro':0.4, 'Py':0.1},     # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
           'MBC_3':{'Fast':0.5,'Slow':0.3,'Necro':0.6, 'Py':0.1},
           'MBC_4':{'Fast':0.5,'Slow':0.3,'Necro':0.6, 'Py':0.1}},
    'tProtected':75.0,        # Protected C turnover time (years). This is the time scale for which protected C is released back to unprotected state.
    'protection_rate':{'Fast':0.0,'Slow':0.001,'Necro':0, 'Py':0}, # Protected carbon formation rate (year-1). Higher number means more will become protected. Can be modified as a function of soil texture/mineralogy to represent different sorption potentials
    'new_resp_units':True,   # At some point I changed the units of vmaxref to be normalized for other factors so they are actually in year-1 units. Leave this as True values that are easier to interpret.
}

envir_vals={'thetamin': array(0.5),
              'thetamax': array(0.7),
              'porosity': array(0.4)}

import Microbial_CORPSE_array
Microbial_CORPSE_array.check_params(params)


# This section is setting up different initial values and parameters for different simulations representing microbial community traits
# Set up empty python dictionaries
initvals={}
paramsets={}
envir_params={}

 
# Relative amount of total C in Sandy vs. Org soils
#    Sandy = 5 g - based on total C (g) estimates from 2022 work
#    Org = 15 g 

# Relative amount of C in Fast vs. Slow C pools
#    Fast = 10%; Slow = 90% of total C based on 2 pool C model (2022)

# Relative amount of PyC
#    5-10% of total C in boreal forest soils based on Soucemarianadin et al., 2014, Geoderma

# Division of C between Necromass and fast C
#     ?

totalC_sandy_unburned = 6
totalC_org_unburned = 15
totalC_sandy_burned = 4.75
totalC_org_burned=14.8

# Initial pools for soil organic matter simulation
# There are three kinds of chemically-defined C (Fast, slow, and microbial necromass). "Fast" has higher maximum decomposition rate and microbial CUE
# Each C type can be in a protected or unprotected state. When protected, it is not subject to microbial decomposition
SOM_init={'CO2': array(0.0),    # Cumulative C-CO2 from microbial respiration
 # To set number of microbial pools, modify "microbial_pools" list in Microbial_CORPSE_array.py
 'MBC_1': array(totalC_sandy_unburned*0.1), # Active, living microbial biomass, 
 'MBC_2': array(totalC_sandy_unburned*0.0),  # MBC_2 = 1% of MBC_1
 'MBC_3': array(0.00),
 'MBC_4': array(0.00),
 'pFastC': array(0.0),         # Protected fast-decomposing C
 'pNecroC': array(0.0),        # Protected microbial necromass C
 'pSlowC': array(0.0),         # Protected slow-decomposing C
 'pPyC': array(0.0),            # Protected PyC
 'uFastC': array(totalC_sandy_unburned*0.089),          # Unprotected fast-decomposing C, 
 'uNecroC': array(totalC_sandy_unburned*0.001),        # Unprotected microbial necromass C, 
 'uSlowC': array(totalC_sandy_unburned*0.85),         # Unprotected slow-decomposing C, 
 'uPyC': array(totalC_sandy_unburned*0.05)}            # Unprotected PyC

# Microbial community, no burn
initvals['no burn sandy soil']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
paramsets['no burn sandy soil']=copy.deepcopy(params)
envir_params['no burn sandy soil']=copy.deepcopy(envir_vals)

# Sandy soil, high severity burn
initvals['high sev burn sandy soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn sandy soil']['MBC_1']=array(totalC_sandy_burned*0.02) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn sandy soil']['MBC_2']=array(totalC_sandy_burned*0.0)
initvals['high sev burn sandy soil']['MBC_3']=array(0.0)
initvals['high sev burn sandy soil']['MBC_4']=array(0.0)
initvals['high sev burn sandy soil']['uFastC']=array(totalC_sandy_burned*0.045) 
initvals['high sev burn sandy soil']['uNecroC']=array(totalC_sandy_burned*0.0031) 
initvals['high sev burn sandy soil']['uSlowC']=array(totalC_sandy_burned*0.86) # Changed value from 80 to 100 based on calibration results
initvals['high sev burn sandy soil']['uPyC']=array(totalC_sandy_burned*0.1)
paramsets['high sev burn sandy soil']=copy.deepcopy(params)

envir_params['high sev burn sandy soil']=copy.deepcopy(envir_vals)
envir_params['high sev burn sandy soil']['thetamin']=array(0.5)
envir_params['high sev burn sandy soil']['thetamax']=array(0.7)
envir_params['high sev burn sandy soil']['porosity']=array(0.4)


# Org. soil, no burn
initvals['no burn org soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['no burn org soil']['MBC_1']=array(totalC_org_unburned*0.1) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['no burn org soil']['MBC_2']=array(totalC_org_unburned*0.0)
initvals['no burn org soil']['MBC_3']=array(0.0)
initvals['no burn org soil']['MBC_4']=array(0.0)
initvals['no burn org soil']['uFastC']=array(totalC_org_unburned*0.10) 
initvals['no burn org soil']['uNecroC']=array(totalC_org_unburned*0.001) 
initvals['no burn org soil']['uSlowC']=array(totalC_org_unburned*0.80) # Changed value from 86.0 to 60 based on calibration results
initvals['no burn org soil']['uPyC']=array(totalC_org_unburned*0.05)
paramsets['no burn org soil']=copy.deepcopy(params)
envir_params['no burn org soil']=copy.deepcopy(envir_params)
envir_params['no burn org soil']['thetamin']=array(0.5)
envir_params['no burn org soil']['thetamax']=array(0.7)
envir_params['no burn org soil']['porosity']=array(0.9)



# Org soil, high severity burn
initvals['high sev burn org soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn org soil']['MBC_1']=array(totalC_org_burned*0.02) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn org soil']['MBC_2']=array(totalC_org_burned*0.0)
initvals['high sev burn org soil']['MBC_3']=array(0.0)
initvals['high sev burn org soil']['MBC_4']=array(0.0)
initvals['high sev burn org soil']['uFastC']=array(totalC_org_burned*0.05) 
initvals['high sev burn org soil']['uNecroC']=array(totalC_org_burned*0.002) 
initvals['high sev burn org soil']['uSlowC']=array(totalC_org_burned*0.84) # C
initvals['high sev burn org soil']['uPyC']=array(totalC_org_burned*0.1)
paramsets['high sev burn org soil']=copy.deepcopy(params)
envir_params['high sev burn org soil']=copy.deepcopy(envir_params)
envir_params['high sev burn org soil']['thetamin']=array(0.5)
envir_params['high sev burn org soil']['thetamax']=array(0.7)
envir_params['high sev burn org soil']['porosity']=array(0.4)


from numpy import where
num_micro_pools=0
num_micro_pools=where(SOM_init['MBC_1']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_2']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_3']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_4']>0, num_micro_pools+1, num_micro_pools+0)   

### Should also include a function that makes sure that each simulation is using the same number of microbial pools!

# Set up a data structure to hold the results of the different simulations
results={}
# Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
# Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
# Inputs are empty because this is running as an incubation without any constant inputs of C



for functype in initvals:
    results[functype] = Microbial_CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=envir_params[functype]['thetamin'],
                                                      thetamax=envir_params[functype]['thetamax'],
                                            times=t,inputs={},clay=2.5,initvals=initvals[functype],params=paramsets[functype])


