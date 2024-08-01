import pandas as pd
import CORPSE_solvers
from numpy import array,arange
import copy


# Initial pools for soil organic matter simulation
# There are three kinds of chemically-defined C (Fast, slow, and microbial necromass). "Fast" has higher maximum decomposition rate and microbial CUE
# Each C type can be in a protected or unprotected state. When protected, it is not subject to microbial decomposition
SOM_init={'CO2': array(0.0),    # Cumulative C-CO2 from microbial respiration
 'livingMicrobeC': array(5.0), # Active, living microbial biomass
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
    'vmaxref':{'Fast':8.0,'Slow':0.2,'Necro':4.0, 'Py':0.1}, #  Relative maximum enzymatic decomp rates for each C type (year-1)
    'Ea':{'Fast':5e3,'Slow':30e3,'Necro':5e3, 'Py':35e3},      # Activation energy (controls T dependence)
    'kC':{'Fast':0.01,'Slow':0.01,'Necro':0.01, 'Py':0.01},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    'gas_diffusion_exp':0.6,  # Determines suppression of decomposition at high soil moisture
    'substrate_diffusion_exp':1.5,   # Controls suppression of decomp at low soil moisture
    'minMicrobeC':1e-3,       # Minimum microbial biomass (fraction of total C). Prevents microbes from going extinct and allows some slow decomposition under adverse conditions
    'Tmic':0.25,              # Microbial biomass mean lifetime (years)
    'et':0.6,                 # Fraction of microbial biomass turnover (death) that goes to necromass instead of being immediately mineralized to CO2
    'eup':{'Fast':0.8,'Slow':0.4,'Necro':0.8, 'Py':0.1},     # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
    'tProtected':75.0,        # Protected C turnover time (years). This is the time scale for which protected C is released back to unprotected state.
    'protection_rate':{'Fast':0.0,'Slow':0.001,'Necro':0, 'Py':0}, # Protected carbon formation rate (year-1). Higher number means more will become protected. Can be modified as a function of soil texture/mineralogy to represent different sorption potentials
    'new_resp_units':True,   # At some point I changed the units of vmaxref to be normalized for other factors so they are actually in year-1 units. Leave this as True values that are easier to interpret.
}


envir_vals={'thetamin': array(0.5),
              'thetamax': array(0.7),
              'porosity': array(0.4)}



import CORPSE_array
CORPSE_array.check_params(params)

# This makes an array of all the model time steps (in units of years). In this case, it starts at zero, ends at 120 days, and has a time step of one day
t=arange(0,70/365,1/365)

# This section is setting up different initial values and parameters for different simulations representing microbial community traits
# Here we set up an empty python dictionary to hold the different sets of parameters and initial values
initvals={}
paramsets={}
envir_params={}

# Microbial community, no burn
initvals['no burn sandy soil']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
paramsets['no burn sandy soil']=copy.deepcopy(params)
envir_params['no burn sandy soil']=copy.deepcopy(envir_vals)

# Sandy soil, high severity burn
initvals['high sev burn sandy soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn sandy soil']['livingMicrobeC']=array(0.6) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn sandy soil']['uFastC']=array(0.75) 
initvals['high sev burn sandy soil']['uNecroC']=array(0.5) 
initvals['high sev burn sandy soil']['uSlowC']=array(80) # Changed value from 80 to 100 based on calibration results
initvals['high sev burn sandy soil']['uPyC']=array(15)
paramsets['high sev burn sandy soil']=copy.deepcopy(params)
paramsets['high sev burn sandy soil']['vmaxref']['Fast']=75.0     # 
paramsets['high sev burn sandy soil']['vmaxref']['Slow']=0.15
paramsets['high sev burn sandy soil']['vmaxref']['Necro']=75.0
paramsets['high sev burn sandy soil']['vmaxref']['Py']=0.05
paramsets['high sev burn sandy soil']['eup']['Fast'] = 0.5        # 
paramsets['high sev burn sandy soil']['eup']['Slow'] = 0.3
paramsets['high sev burn sandy soil']['eup']['Necro'] = 0.5
paramsets['high sev burn sandy soil']['eup']['Py'] = 0.1
envir_params['high sev burn sandy soil']=copy.deepcopy(envir_vals)
envir_params['high sev burn sandy soil']['thetamin']=array(0.5)
envir_params['high sev burn sandy soil']['thetamax']=array(0.7)
envir_params['high sev burn sandy soil']['porosity']=array(0.4)

# Org. soil, no burn
initvals['no burn org soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['no burn org soil']['livingMicrobeC']=array(6.0) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['no burn org soil']['uFastC']=array(6.0) 
initvals['no burn org soil']['uNecroC']=array(1.2) 
initvals['no burn org soil']['uSlowC']=array(72) # Changed value from 86.0 to 60 based on calibration results
initvals['no burn org soil']['uPyC']=array(5.0)
paramsets['no burn org soil']=copy.deepcopy(params)
paramsets['no burn org soil']['vmaxref']['Fast']=7.25     # 
paramsets['no burn org soil']['vmaxref']['Slow']=0.1
paramsets['no burn org soil']['vmaxref']['Necro']=7.25
paramsets['no burn org soil']['vmaxref']['Py']=0.1
paramsets['no burn org soil']['eup']['Fast'] = 0.9        # 
paramsets['no burn org soil']['eup']['Slow'] = 0.4
paramsets['no burn org soil']['eup']['Necro'] = 0.7
paramsets['no burn org soil']['eup']['Py'] = 0.1
envir_params['no burn org soil']=copy.deepcopy(envir_params)
envir_params['no burn org soil']['thetamin']=array(0.5)
envir_params['no burn org soil']['thetamax']=array(0.7)
envir_params['no burn org soil']['porosity']=array(0.9)


# Org soil, high severity burn
initvals['high sev burn org soil']=copy.deepcopy(SOM_init)       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn org soil']['livingMicrobeC']=array(0.6) # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn org soil']['uFastC']=array(0.6) 
initvals['high sev burn org soil']['uNecroC']=array(0.6) 
initvals['high sev burn org soil']['uSlowC']=array(86) # Changed value from 86 to 60 based on calibration results
initvals['high sev burn org soil']['uPyC']=array(12)
paramsets['high sev burn org soil']=copy.deepcopy(params)
paramsets['high sev burn org soil']['vmaxref']['Fast']=75.0     # 
paramsets['high sev burn org soil']['vmaxref']['Slow']=0.15
paramsets['high sev burn org soil']['vmaxref']['Necro']=75.0
paramsets['high sev burn org soil']['vmaxref']['Py']=0.1
paramsets['high sev burn org soil']['eup']['Fast'] = 0.4        # 
paramsets['high sev burn org soil']['eup']['Slow'] = 0.3
paramsets['high sev burn org soil']['eup']['Necro'] = 0.4
paramsets['high sev burn org soil']['eup']['Py'] = 0.1
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
    results[functype] = CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=envir_params[functype]['thetamin'],
                                                      thetamax=envir_params[functype]['thetamax'],
                                            times=t,inputs={},clay=2.5,initvals=initvals[functype],params=paramsets[functype])




## microbial pools:
# # Resistant: Low biomass loss, fast growth
# # Fast-growing survivor
# initvals['Fast-growing survivor']=copy.deepcopy(SOM_init)      # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
# initvals['Fast-growing survivor']['livingMicrobeC']=array(0.01) # Start with a high amount of initial microbial biomass. Assumes this simulation starts right after the fire
# paramsets['Fast-growing survivor']=copy.deepcopy(params)
# paramsets['Fast-growing survivor']['vmaxref']['Fast']=20.0     # Higher potential decomposition rate for more labile C
# paramsets['Fast-growing survivor']['eup']['Fast'] = 0.2
# paramsets['Fast-growing survivor']['eup']['Necro'] = 0.2


# # Susceptible: High biomass loss, slow growth
# # Fire susceptible
# initvals['Fire susceptible']=copy.deepcopy(SOM_init)
# initvals['Fire susceptible']['livingMicrobeC']=array(0.001)    # Low initial microbial biomass
# paramsets['Fire susceptible']=copy.deepcopy(params)
# paramsets['Fire susceptible']['vmaxref']['Fast']=0.2          # Slower maximum decomposition rate for labile C
# paramsets['Fire susceptible']['eup']['Fast'] = 0.5
# paramsets['Fire susceptible']['eup']['Necro'] = 0.5


# # Recovering: High biomass loss, fast growth
# # Post-fire rebounder
# initvals['Post-fire rebounder']=copy.deepcopy(SOM_init)
# initvals['Post-fire rebounder']['livingMicrobeC']=array(0.001) # Low initial microbial biomass
# paramsets['Post-fire rebounder']=copy.deepcopy(params)
# paramsets['Post-fire rebounder']['vmaxref']['Fast']=20.0      # Very fast potential decomposition/growth rate
# paramsets['Post-fire rebounder']['eup']['Fast'] = 0.2
# paramsets['Post-fire rebounder']['eup']['Necro'] = 0.2


# # Resilient: Low biomass loss, slow growth
# # Slow-growing survivor
# initvals['Slow-growing survivor']=copy.deepcopy(SOM_init)
# initvals['Slow-growing survivor']['livingMicrobeC']=array(0.01) # High initial microbial biomass
# paramsets['Slow-growing survivor']=copy.deepcopy(params)
# paramsets['Slow-growing survivor']['vmaxref']['Fast']=0.2      # Slower decomposition/growth rate
# paramsets['Slow-growing survivor']['eup']['Fast'] = 0.5
# paramsets['Slow-growing survivor']['eup']['Necro'] = 0.5


