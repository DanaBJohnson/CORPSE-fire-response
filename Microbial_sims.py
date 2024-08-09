import pandas as pd
import Microbial_CORPSE_solvers
from numpy import array,arange
import copy


# Initial pools for soil organic matter simulation
# There are three kinds of chemically-defined C (Fast, slow, and microbial necromass). "Fast" has higher maximum decomposition rate and microbial CUE
# Each C type can be in a protected or unprotected state. When protected, it is not subject to microbial decomposition
SOM_init={'CO2': array(0.0),    # Cumulative C-CO2 from microbial respiration
 'MBC_1': array(4.95), # Active, living microbial biomass, 
 'MBC_2': array(0.05), # To set number of microbial pools, modify "microbial_pools" list in Microbial_CORPSE_array.py
 'MBC_3': array(0.00),
 'MBC_4': array(0.00),
 'pFastC': array(0.0),         # Protected fast-decomposing C
 'pNecroC': array(0.0),        # Protected microbial necromass C
 'pSlowC': array(0.0),         # Protected slow-decomposing C
 'pPyC': array(0.0),            # Protected PyC
 'uFastC': array(3.0),          # Unprotected fast-decomposing C, Changed value from 10 to 11 based on calibration results
 'uNecroC': array(3.0),        # Unprotected microbial necromass C, 
 'uSlowC': array(88),         # Unprotected slow-decomposing C, Changed value from 86 to 60.2 based on calibration results
 'uPyC': array(4.0)}            # Unprotected PyC
 
# Parameters controlling the model (sandy soil, no burn)
params={
    #'vmaxref':{'Fast':8.0,'Slow':0.2,'Necro':4.0, 'Py':0.1}, #  Relative maximum enzymatic decomp rates for each C type (year-1)
    'vmaxref_1':{'Fast':8.0,'Slow':0.2,'Necro':4.0, 'Py':0.1}, # Slower-growing microbial pool should have lower Vmax, Button et al. 1991
    'vmaxref_2':{'Fast':25.0,'Slow':0.02,'Necro':12.0, 'Py':0.01}, # Faster-growing microbial pool should have higher Vmax
    'vmaxref_3':{'Fast':0.0,'Slow':0.0,'Necro':0.0, 'Py':0.0}, # 
    'vmaxref_4':{'Fast':0.0,'Slow':0.0,'Necro':0.0, 'Py':0.0}, # 
    'Ea':{'Fast':5e3,'Slow':30e3,'Necro':5e3, 'Py':35e3},      # Activation energy (controls T dependence)
    'kC_1':{'Fast':0.01,'Slow':0.01,'Necro':0.01, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    'kC_2':{'Fast':0.01,'Slow':0.01,'Necro':0.01, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    'kC_3':{'Fast':0.1,'Slow':0.01,'Necro':0.1, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    'kC_4':{'Fast':0.1,'Slow':0.01,'Necro':0.1, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    #### Should kC vary for microbial pools??
    'gas_diffusion_exp':0.6,  # Determines suppression of decomposition at high soil moisture
    'substrate_diffusion_exp':1.5,   # Controls suppression of decomp at low soil moisture
    'minMicrobeC':1e-4,       # Minimum microbial biomass (fraction of total C). Prevents microbes from going extinct and allows some slow decomposition under adverse conditions
    'Tmic':0.25,              # Microbial biomass mean lifetime (years)
    'et_1':0.6,                 # Fraction of microbial biomass turnover (death) that goes to necromass instead of being immediately mineralized to CO2
    'et_2':0.6,                 # Fraction of microbial biomass turnover (death) that goes to necromass instead of being immediately mineralized to CO2
    'et_3':0.6, 
    'et_4':0.6, 
    'eup_1':{'Fast':0.6,'Slow':0.4,'Necro':0.6, 'Py':0.2},     # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
    'eup_2':{'Fast':0.7,'Slow':0.3,'Necro':0.7, 'Py':0.1},     # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
    'eup_3':{'Fast':0.5,'Slow':0.3,'Necro':0.6, 'Py':0.1},
    'eup_4':{'Fast':0.5,'Slow':0.3,'Necro':0.6, 'Py':0.1},
    'tProtected':75.0,        # Protected C turnover time (years). This is the time scale for which protected C is released back to unprotected state.
    'protection_rate':{'Fast':0.0,'Slow':0.001,'Necro':0, 'Py':0}, # Protected carbon formation rate (year-1). Higher number means more will become protected. Can be modified as a function of soil texture/mineralogy to represent different sorption potentials
    'new_resp_units':True,   # At some point I changed the units of vmaxref to be normalized for other factors so they are actually in year-1 units. Leave this as True values that are easier to interpret.
}

envir_vals={'thetamin': array(0.5),
              'thetamax': array(0.7),
              'porosity': array(0.4)}

from numpy import where
num_micro_pools=0
num_micro_pools=where(SOM_init['MBC_1']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_2']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_3']>0, num_micro_pools+1, num_micro_pools+0)   
num_micro_pools=where(SOM_init['MBC_4']>0, num_micro_pools+1, num_micro_pools+0)   

### Should also include a function that makes sure that each simulation is using the same number of microbial pools!

import Microbial_CORPSE_array
Microbial_CORPSE_array.check_params(params)

# This makes an array of all the model time steps (in units of years). In this case, it starts at zero, ends at 120 days, and has a time step of one day
t=arange(0,70/365,1/365)

# This section is setting up different initial values and parameters for different simulations representing microbial community traits
# Set up empty python dictionaries
initvals={}
paramsets={}
envir_params={}

# Microbial community, no burn
initvals['no burn sandy soil']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
paramsets['no burn sandy soil']=copy.deepcopy(params)
envir_params['no burn sandy soil']=copy.deepcopy(envir_vals)

# Sandy soil, high severity burn
initvals['high sev burn sandy soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn sandy soil']['MBC_1']=array(0.44) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn sandy soil']['MBC_2']=array(0.15)
initvals['high sev burn sandy soil']['MBC_3']=array(0.0)
initvals['high sev burn sandy soil']['MBC_4']=array(0.0)
initvals['high sev burn sandy soil']['uFastC']=array(0.75) 
initvals['high sev burn sandy soil']['uNecroC']=array(1.0) 
initvals['high sev burn sandy soil']['uSlowC']=array(80) # Changed value from 80 to 100 based on calibration results
initvals['high sev burn sandy soil']['uPyC']=array(15)
paramsets['high sev burn sandy soil']=copy.deepcopy(params)
#paramsets['high sev burn sandy soil']['vmaxref_1']['Fast']=8.0 # MBC 1, vmax
#paramsets['high sev burn sandy soil']['vmaxref_1']['Slow']=0.2
#paramsets['high sev burn sandy soil']['vmaxref_1']['Necro']=4.0
#paramsets['high sev burn sandy soil']['vmaxref_1']['Py']=0.1
paramsets['high sev burn sandy soil']['vmaxref_2']['Fast']=paramsets['high sev burn sandy soil']['vmaxref_1']['Fast']*10   # 
paramsets['high sev burn sandy soil']['vmaxref_2']['Slow']=paramsets['high sev burn sandy soil']['vmaxref_1']['Slow']*0.1
paramsets['high sev burn sandy soil']['vmaxref_2']['Necro']=paramsets['high sev burn sandy soil']['vmaxref_1']['Necro']*10
paramsets['high sev burn sandy soil']['vmaxref_2']['Py']=paramsets['high sev burn sandy soil']['vmaxref_1']['Py']*0.1
paramsets['high sev burn sandy soil']['vmaxref_3']['Fast']=75.0 # MBC 3, vmax
paramsets['high sev burn sandy soil']['vmaxref_3']['Slow']=0.15
paramsets['high sev burn sandy soil']['vmaxref_3']['Necro']=75.0
paramsets['high sev burn sandy soil']['vmaxref_3']['Py']=0.05
paramsets['high sev burn sandy soil']['vmaxref_4']['Fast']=75.0 # MBC 4, vmax
paramsets['high sev burn sandy soil']['vmaxref_4']['Slow']=0.15
paramsets['high sev burn sandy soil']['vmaxref_4']['Necro']=75.0
paramsets['high sev burn sandy soil']['vmaxref_4']['Py']=0.05
paramsets['high sev burn sandy soil']['kC_1']['Fast']=0.01 # MBC 1, kC
paramsets['high sev burn sandy soil']['kC_2']['Slow']=0.01 # MBC 2, kC
paramsets['high sev burn sandy soil']['kC_3']['Necro']=0.01 # MBC 3, kC
paramsets['high sev burn sandy soil']['kC_4']['Py']=0.01 # MBC 4, kC
#paramsets['high sev burn sandy soil']['eup_1']['Fast'] = 0.5 # MBC 1, CUE
#paramsets['high sev burn sandy soil']['eup_1']['Slow'] = 0.3
#paramsets['high sev burn sandy soil']['eup_1']['Necro'] = 0.5
#paramsets['high sev burn sandy soil']['eup_1']['Py'] = 0.05
paramsets['high sev burn sandy soil']['eup_2']['Fast'] = paramsets['high sev burn sandy soil']['eup_1']['Fast']*1.5 # MBC 2, CUE
paramsets['high sev burn sandy soil']['eup_2']['Slow'] = paramsets['high sev burn sandy soil']['eup_1']['Slow']
paramsets['high sev burn sandy soil']['eup_2']['Necro'] = paramsets['high sev burn sandy soil']['eup_1']['Necro']*1.5
paramsets['high sev burn sandy soil']['eup_2']['Py'] = paramsets['high sev burn sandy soil']['eup_1']['Py']
paramsets['high sev burn sandy soil']['eup_3']['Fast'] = 0.4 # MBC 3, CUE
paramsets['high sev burn sandy soil']['eup_3']['Slow'] = 0.2
paramsets['high sev burn sandy soil']['eup_3']['Necro'] = 0.4
paramsets['high sev burn sandy soil']['eup_3']['Py'] = 0.1
paramsets['high sev burn sandy soil']['eup_4']['Fast'] = 0.4 # MBC 4, CUE
paramsets['high sev burn sandy soil']['eup_4']['Slow'] = 0.2
paramsets['high sev burn sandy soil']['eup_4']['Necro'] = 0.4
paramsets['high sev burn sandy soil']['eup_4']['Py'] = 0.1
envir_params['high sev burn sandy soil']=copy.deepcopy(envir_vals)
envir_params['high sev burn sandy soil']['thetamin']=array(0.5)
envir_params['high sev burn sandy soil']['thetamax']=array(0.7)
envir_params['high sev burn sandy soil']['porosity']=array(0.4)


# Org. soil, no burn
initvals['no burn org soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['no burn org soil']['MBC_1']=array(5.98) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['no burn org soil']['MBC_2']=array(0.02)
initvals['no burn org soil']['MBC_3']=array(0.0)
initvals['no burn org soil']['MBC_4']=array(0.0)
initvals['no burn org soil']['uFastC']=array(6.0) 
initvals['no burn org soil']['uNecroC']=array(1.2) 
initvals['no burn org soil']['uSlowC']=array(72) # Changed value from 86.0 to 60 based on calibration results
initvals['no burn org soil']['uPyC']=array(5.0)
paramsets['no burn org soil']=copy.deepcopy(params)
paramsets['no burn org soil']['vmaxref_1']['Fast']=2.5    # 
paramsets['no burn org soil']['vmaxref_1']['Slow']=0.1
paramsets['no burn org soil']['vmaxref_1']['Necro']=2.5
paramsets['no burn org soil']['vmaxref_1']['Py']=0.1
paramsets['no burn org soil']['vmaxref_2']['Fast']=paramsets['no burn org soil']['vmaxref_1']['Fast']*10 # 
paramsets['no burn org soil']['vmaxref_2']['Slow']=paramsets['no burn org soil']['vmaxref_1']['Slow']*0.1
paramsets['no burn org soil']['vmaxref_2']['Necro']=paramsets['no burn org soil']['vmaxref_1']['Necro']*10
paramsets['no burn org soil']['vmaxref_2']['Py']=paramsets['no burn org soil']['vmaxref_1']['Py']*0.1
paramsets['no burn org soil']['vmaxref_3']['Fast']=7.25     # 
paramsets['no burn org soil']['vmaxref_3']['Slow']=0.1
paramsets['no burn org soil']['vmaxref_3']['Necro']=7.25
paramsets['no burn org soil']['vmaxref_3']['Py']=0.1
paramsets['no burn org soil']['vmaxref_4']['Fast']=7.25     # 
paramsets['no burn org soil']['vmaxref_4']['Slow']=0.1
paramsets['no burn org soil']['vmaxref_4']['Necro']=7.25
paramsets['no burn org soil']['vmaxref_4']['Py']=0.1
paramsets['no burn org soil']['kC_1']['Fast']=0.01     # 
paramsets['no burn org soil']['kC_2']['Slow']=0.01
paramsets['no burn org soil']['kC_3']['Necro']=0.01
paramsets['no burn org soil']['kC_4']['Py']=0.01
paramsets['no burn org soil']['eup_1']['Fast'] = 0.5        # 
paramsets['no burn org soil']['eup_1']['Slow'] = 0.3
paramsets['no burn org soil']['eup_1']['Necro'] = 0.5
paramsets['no burn org soil']['eup_1']['Py'] = 0.1
paramsets['no burn org soil']['eup_2']['Fast'] = paramsets['no burn org soil']['eup_1']['Fast']*1.5
paramsets['no burn org soil']['eup_2']['Slow'] = paramsets['no burn org soil']['eup_1']['Slow']
paramsets['no burn org soil']['eup_2']['Necro'] = paramsets['no burn org soil']['eup_1']['Necro']*1.5
paramsets['no burn org soil']['eup_2']['Py'] = paramsets['no burn org soil']['eup_1']['Py']
paramsets['no burn org soil']['eup_3']['Fast'] = 0.8        # 
paramsets['no burn org soil']['eup_3']['Slow'] = 0.3
paramsets['no burn org soil']['eup_3']['Necro'] = 0.6
paramsets['no burn org soil']['eup_3']['Py'] = 0.05
paramsets['no burn org soil']['eup_4']['Fast'] = 0.9        # 
paramsets['no burn org soil']['eup_4']['Slow'] = 0.4
paramsets['no burn org soil']['eup_4']['Necro'] = 0.7
paramsets['no burn org soil']['eup_4']['Py'] = 0.1
envir_params['no burn org soil']=copy.deepcopy(envir_params)
envir_params['no burn org soil']['thetamin']=array(0.5)
envir_params['no burn org soil']['thetamax']=array(0.7)
envir_params['no burn org soil']['porosity']=array(0.9)


# Org soil, high severity burn
initvals['high sev burn org soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn org soil']['MBC_1']=array(0.51) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn org soil']['MBC_2']=array(0.09)
initvals['high sev burn org soil']['MBC_3']=array(0.0)
initvals['high sev burn org soil']['MBC_4']=array(0.0)
initvals['high sev burn org soil']['uFastC']=array(1.) 
initvals['high sev burn org soil']['uNecroC']=array(1.0) 
initvals['high sev burn org soil']['uSlowC']=array(70) # Changed value from 86 to 60 based on calibration results
initvals['high sev burn org soil']['uPyC']=array(6)
paramsets['high sev burn org soil']=copy.deepcopy(params)
paramsets['high sev burn org soil']['vmaxref_1']['Fast']=paramsets['no burn org soil']['vmaxref_1']['Fast'] 
paramsets['high sev burn org soil']['vmaxref_1']['Slow']=paramsets['no burn org soil']['vmaxref_1']['Slow']
paramsets['high sev burn org soil']['vmaxref_1']['Necro']=paramsets['no burn org soil']['vmaxref_1']['Necro']
paramsets['high sev burn org soil']['vmaxref_1']['Py']=paramsets['no burn org soil']['vmaxref_1']['Py']
paramsets['high sev burn org soil']['vmaxref_2']['Fast']=paramsets['high sev burn org soil']['vmaxref_1']['Fast']*10 # 
paramsets['high sev burn org soil']['vmaxref_2']['Slow']=paramsets['high sev burn org soil']['vmaxref_1']['Slow']*0.1
paramsets['high sev burn org soil']['vmaxref_2']['Necro']=paramsets['high sev burn org soil']['vmaxref_1']['Necro']*10
paramsets['high sev burn org soil']['vmaxref_2']['Py']=paramsets['high sev burn org soil']['vmaxref_1']['Py']*0.1
paramsets['high sev burn org soil']['vmaxref_3']['Fast']=75.0     # 
paramsets['high sev burn org soil']['vmaxref_3']['Slow']=0.15
paramsets['high sev burn org soil']['vmaxref_3']['Necro']=75.0
paramsets['high sev burn org soil']['vmaxref_3']['Py']=0.1
paramsets['high sev burn org soil']['vmaxref_4']['Fast']=75.0     # 
paramsets['high sev burn org soil']['vmaxref_4']['Slow']=0.15
paramsets['high sev burn org soil']['vmaxref_4']['Necro']=75.0
paramsets['high sev burn org soil']['vmaxref_4']['Py']=0.1
paramsets['high sev burn org soil']['kC_1']['Fast']=0.01     # 
paramsets['high sev burn org soil']['kC_2']['Slow']=0.01
paramsets['high sev burn org soil']['kC_3']['Necro']=0.01
paramsets['high sev burn org soil']['kC_4']['Py']=0.01
paramsets['high sev burn org soil']['eup_1']['Fast'] = paramsets['no burn org soil']['eup_1']['Fast']
paramsets['high sev burn org soil']['eup_1']['Slow'] = paramsets['no burn org soil']['eup_1']['Slow']
paramsets['high sev burn org soil']['eup_1']['Necro'] = paramsets['no burn org soil']['eup_1']['Necro']
paramsets['high sev burn org soil']['eup_1']['Py'] = paramsets['no burn org soil']['eup_1']['Py']
paramsets['high sev burn org soil']['eup_2']['Fast'] = paramsets['high sev burn org soil']['eup_1']['Fast']*1.5
paramsets['high sev burn org soil']['eup_2']['Slow'] = paramsets['high sev burn org soil']['eup_1']['Slow']
paramsets['high sev burn org soil']['eup_2']['Necro'] = paramsets['high sev burn org soil']['eup_1']['Necro']*1.5
paramsets['high sev burn org soil']['eup_2']['Py'] = paramsets['high sev burn org soil']['eup_1']['Py']
paramsets['high sev burn org soil']['eup_3']['Fast'] = 0.3        # 
paramsets['high sev burn org soil']['eup_3']['Slow'] = 0.2
paramsets['high sev burn org soil']['eup_3']['Necro'] = 0.3
paramsets['high sev burn org soil']['eup_3']['Py'] = 0.05
paramsets['high sev burn org soil']['eup_4']['Fast'] = 0.4        # 
paramsets['high sev burn org soil']['eup_4']['Slow'] = 0.3
paramsets['high sev burn org soil']['eup_4']['Necro'] = 0.4
paramsets['high sev burn org soil']['eup_4']['Py'] = 0.1
envir_params['high sev burn org soil']=copy.deepcopy(envir_params)
envir_params['high sev burn org soil']['thetamin']=array(0.5)
envir_params['high sev burn org soil']['thetamax']=array(0.7)
envir_params['high sev burn org soil']['porosity']=array(0.9)

# Set up a data structure to hold the results of the different simulations
results={}
# Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
# Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
# Inputs are empty because this is running as an incubation without any constant inputs of C



for functype in initvals:
    results[functype] = Microbial_CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=envir_params[functype]['thetamin'],
                                                      thetamax=envir_params[functype]['thetamax'],
                                            times=t,inputs={},clay=2.5,initvals=initvals[functype],params=paramsets[functype])


