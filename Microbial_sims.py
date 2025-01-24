import pandas as pd
from numpy import where
import CORPSE_array
import CORPSE_solvers
from numpy import array, arange
import copy

# This makes an array of all the model time steps (in units of years). In this case, it starts at zero, ends at 120 days, and has a time step of one day
t = arange(0, 365/365, 1/365)

# Parameters controlling the model
params = {
    # 'vmaxref':{'Fast':8.0,'Slow':0.2,'Necro':4.0, 'Py':0.1}, #  Relative maximum enzymatic decomp rates for each C type (year-1)
    'vmaxref': {'MBC_1': {'Fast': 6.9, 'Slow': 0.11, 'Necro': 7.0, 'Py': 0.1},
                'MBC_2': {'Fast': 19.2, 'Slow': 0.0064, 'Necro': 45.0, 'Py': 0.01},
                'MBC_3': {'Fast': 0.0, 'Slow': 0.0, 'Necro': 0.0, 'Py': 0.0},
                'MBC_4': {'Fast': 0.0, 'Slow': 0.0, 'Necro': 0.0, 'Py': 0.0}},
    # Activation energy (controls T dependence)
    'Ea': {'Fast': 30e3, 'Slow': 30e3, 'Necro': 30e3, 'Py': 35e3},
    'kC': {'MBC_1': {'Fast': 0.0085, 'Slow': 0.02, 'Necro': 0.01, 'Py': 0.02},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
           # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
           'MBC_2': {'Fast': 0.01, 'Slow': 0.04, 'Necro': 0.01, 'Py': 0.04},
           # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
           'MBC_3': {'Fast': 0.1, 'Slow': 0.01, 'Necro': 0.1, 'Py': 0.01},
           'MBC_4': {'Fast': 0.1, 'Slow': 0.01, 'Necro': 0.1, 'Py': 0.01}},    # Michaelis-Menton half saturation parameter (g microbial biomass/g substrate)
    #       Smaller kC results in more CO2 produced for a given amount of C available to decompose?
    #       Copiotrophs should have larger kC than oligotrophs. Oligotrophs
    #       have higher substrate affinity.
    # Determines suppression of decomposition at high soil moisture
    'gas_diffusion_exp': 0.6,
    # Controls suppression of decomp at low soil moisture
    'substrate_diffusion_exp': 1.5,
    # Minimum microbial biomass (fraction of total C). Prevents microbes from going extinct and allows some slow decomposition under adverse conditions
    'minMicrobeC': {'MBC_1': 1e-3, 'MBC_2': 1e-5, 'MBC_3': 1e-3, 'MBC_4': 1e-3},
    # Microbial biomass mean lifetime (years)
    'Tmic': {'MBC_1': 0.5, 'MBC_2': 0.15, 'MBC_3': 0.25, 'MBC_4': 0.25},
    #       Give MBC 1 a mean lifetime of 6 hours (0.25) and MBC_2 a mean lifetime of 60 minutes.
    # Fraction of microbial biomass turnover (death) that goes to necromass instead of being immediately mineralized to CO2
    'et': {'MBC_1': 0.8, 'MBC_2': 0.8, 'MBC_3': 0.6, 'MBC_4': 0.6}, # Wouldn't we expect all the MBC to become necromass and then decomposition of necromass --> CO2??
    # Microbial carbon use efficiency for each substrate type (fast, slow, necromass). This amount is converted to biomass during decomposition, and the remainder is immediately respired as CO2
    'eup': {'MBC_1': {'Fast': 0.4, 'Slow': 0.3, 'Necro': 0.55, 'Py': 0.15},     # Set based on new CUE data
            'MBC_2': {'Fast': 0.36, 'Slow': 0.1, 'Necro': 0.3, 'Py': 0.05},
            'MBC_3': {'Fast': 0.5, 'Slow': 0.3, 'Necro': 0.6, 'Py': 0.1},
            'MBC_4': {'Fast': 0.5, 'Slow': 0.3, 'Necro': 0.6, 'Py': 0.1}},
    # Protected C turnover time (years). This is the time scale for which protected C is released back to unprotected state.
    'tProtected': 75.0,
    # Protected carbon formation rate (year-1). Higher number means more will become protected. Can be modified as a function of soil texture/mineralogy to represent different sorption potentials
    'protection_rate': {'Fast': 0.0, 'Slow': 0.001, 'Necro': 0, 'Py': 0}, # Rates from previous iteration: 'protection_rate':{'Fast':0.1,'Slow':0.0001,'Necro':1.5}
    # At some point Ben changed the units of vmaxref to be normalized for other factors so they are actually in year-1 units. Leave this as True values that are easier to interpret.
    'new_resp_units': True,
}


envir_vals = {'thetamin': array(0.5),
              'thetamax': array(0.7),
              'porosity': array(0.4)}

CORPSE_array.check_params(params)


# This section is setting up different initial values and parameters for different simulations representing microbial community traits
# Set up empty python dictionaries
initvals = {}
paramsets = {}
envir_params = {}

# Initial total C pool size = from Johnson et al. 2024, Fire Ecology
totalC_sandy_unburned = 5.7
totalC_sandy_low_sev = 5.0
totalC_sandy_high_sev = 4.4

totalC_org_unburned = 14.5
totalC_org_low_sev = 14.7
# average from Johnson et al. 2023 & Johnson et al. 2024 papers
totalC_org_high_sev = 13.45


# Initial pools for soil organic matter simulation
# There are three kinds of chemically-defined C (Fast, slow, and microbial necromass). "Fast" has higher maximum decomposition rate and microbial CUE
# Each C type can be in a protected or unprotected state. When protected, it is not subject to microbial decomposition
SOM_init = {'CO2': array(0.0),    # Cumulative C-CO2 from microbial respiration
            # To set number of microbial pools, modify "microbial_pools" list in CORPSE_array.py
            # Active, living microbial biomass
            'MBC_1': array(totalC_sandy_unburned*0.013*0.99),
            # MBC_2 = 1% of MBC_1
            'MBC_2': array(totalC_sandy_unburned*0.013*0.01),
            'MBC_3': array(0.00),
            'MBC_4': array(0.00),
            'pFastC': array(0.0),         # Protected fast-decomposing C
            'pNecroC': array(0.0),        # Protected microbial necromass C
            'pSlowC': array(0.0),         # Protected slow-decomposing C
            'pPyC': array(0.0),            # Protected PyC
            # Unprotected fast-decomposing C,
            'uFastC': array(totalC_sandy_unburned*0.063),
            # Unprotected microbial necromass C,
            'uNecroC': array(totalC_sandy_unburned*0.018),
            # Unprotected slow-decomposing C,
            'uSlowC': array(totalC_sandy_unburned*0.86),
            'uPyC': array(totalC_sandy_unburned*0.050)}            # Unprotected PyC

# Microbial community, no burn
# Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['no burn Gleysol'] = copy.deepcopy(SOM_init)
paramsets['no burn Gleysol'] = copy.deepcopy(params)
envir_params['no burn Gleysol'] = copy.deepcopy(envir_vals)


# Histosol, no burn
# Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['no burn Histosol'] = copy.deepcopy(SOM_init)
# Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['no burn Histosol']['MBC_1'] = array(totalC_org_unburned*0.0056*0.99)
initvals['no burn Histosol']['MBC_2'] = array(totalC_org_unburned*0.0056*0.01)
initvals['no burn Histosol']['MBC_3'] = array(0.0)
initvals['no burn Histosol']['MBC_4'] = array(0.0)
initvals['no burn Histosol']['uFastC'] = array(totalC_org_unburned*0.077)
initvals['no burn Histosol']['uNecroC'] = array(totalC_org_unburned*0.016)
# Changed value from 86.0 to 60 based on calibration results
initvals['no burn Histosol']['uSlowC'] = array(totalC_org_unburned*0.86)
initvals['no burn Histosol']['uPyC'] = array(totalC_org_unburned*0.045)
paramsets['no burn Histosol'] = copy.deepcopy(params)
envir_params['no burn Histosol'] = copy.deepcopy(envir_params)
envir_params['no burn Histosol']['thetamin'] = array(0.5)
envir_params['no burn Histosol']['thetamax'] = array(0.7)
envir_params['no burn Histosol']['porosity'] = array(0.9)


# Gleysol, low severity burn
# Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['low sev burn Gleysol'] = copy.deepcopy(SOM_init)
# Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['low sev burn Gleysol']['MBC_1'] = array(
    totalC_sandy_low_sev*0.011*0.95)
initvals['low sev burn Gleysol']['MBC_2'] = array(
    totalC_sandy_low_sev*0.011*0.05)
initvals['low sev burn Gleysol']['MBC_3'] = array(0.0)
initvals['low sev burn Gleysol']['MBC_4'] = array(0.0)
initvals['low sev burn Gleysol']['uFastC'] = array(totalC_sandy_low_sev*0.066)
initvals['low sev burn Gleysol']['uNecroC'] = array(totalC_sandy_low_sev*0.015)
# Changed value from 80 to 100 based on calibration results
initvals['low sev burn Gleysol']['uSlowC'] = array(totalC_sandy_low_sev*0.83)
initvals['low sev burn Gleysol']['uPyC'] = array(totalC_sandy_low_sev*0.076)
paramsets['low sev burn Gleysol'] = copy.deepcopy(params)
envir_params['low sev burn Gleysol'] = copy.deepcopy(envir_vals)
envir_params['low sev burn Gleysol']['thetamin'] = array(0.5)
envir_params['low sev burn Gleysol']['thetamax'] = array(0.7)
envir_params['low sev burn Gleysol']['porosity'] = array(0.4)


# Histosol, low severity burn
# Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['low sev burn Histosol'] = copy.deepcopy(SOM_init)
# Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['low sev burn Histosol']['MBC_1'] = array(
    totalC_org_low_sev*0.0050*0.95)
initvals['low sev burn Histosol']['MBC_2'] = array(
    totalC_org_low_sev*0.0050*0.05)
initvals['low sev burn Histosol']['MBC_3'] = array(0.0)
initvals['low sev burn Histosol']['MBC_4'] = array(0.0)
initvals['low sev burn Histosol']['uFastC'] = array(totalC_org_low_sev*0.041)
initvals['low sev burn Histosol']['uNecroC'] = array(totalC_org_low_sev*0.017)
initvals['low sev burn Histosol']['uSlowC'] = array(
    totalC_org_low_sev*0.80)  # C
initvals['low sev burn Histosol']['uPyC'] = array(totalC_org_low_sev*0.14)
paramsets['low sev burn Histosol'] = copy.deepcopy(params)
envir_params['low sev burn Histosol'] = copy.deepcopy(envir_params)
envir_params['low sev burn Histosol']['thetamin'] = array(0.5)
envir_params['low sev burn Histosol']['thetamax'] = array(0.7)
envir_params['low sev burn Histosol']['porosity'] = array(0.4)


# Gleysol, high severity burn
# Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn Gleysol'] = copy.deepcopy(SOM_init)
# Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn Gleysol']['MBC_1'] = array(
    totalC_sandy_high_sev*0.0080*0.85)
initvals['high sev burn Gleysol']['MBC_2'] = array(
    totalC_sandy_high_sev*0.0080*0.15)
initvals['high sev burn Gleysol']['MBC_3'] = array(0.0)
initvals['high sev burn Gleysol']['MBC_4'] = array(0.0)
initvals['high sev burn Gleysol']['uFastC'] = array(
    totalC_sandy_high_sev*0.025)
initvals['high sev burn Gleysol']['uNecroC'] = array(
    totalC_sandy_high_sev*0.016)
# Changed value from 80 to 100 based on calibration results
initvals['high sev burn Gleysol']['uSlowC'] = array(totalC_sandy_high_sev*0.77)
initvals['high sev burn Gleysol']['uPyC'] = array(totalC_sandy_high_sev*0.18)
paramsets['high sev burn Gleysol'] = copy.deepcopy(params)
envir_params['high sev burn Gleysol'] = copy.deepcopy(envir_vals)
envir_params['high sev burn Gleysol']['thetamin'] = array(0.5)
envir_params['high sev burn Gleysol']['thetamax'] = array(0.7)
envir_params['high sev burn Gleysol']['porosity'] = array(0.4)


# Histosol, high severity burn
# Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
initvals['high sev burn Histosol'] = copy.deepcopy(SOM_init)
# Start with low initial microbial biomass. Assumes this simulation starts right after the fire
initvals['high sev burn Histosol']['MBC_1'] = array(
    totalC_org_high_sev*0.0078*0.93)
initvals['high sev burn Histosol']['MBC_2'] = array(
    totalC_org_high_sev*0.0078*0.07)
initvals['high sev burn Histosol']['MBC_3'] = array(0.0)
initvals['high sev burn Histosol']['MBC_4'] = array(0.0)
initvals['high sev burn Histosol']['uFastC'] = array(totalC_org_high_sev*0.018)
initvals['high sev burn Histosol']['uNecroC'] = array(
    totalC_org_high_sev*0.019)
initvals['high sev burn Histosol']['uSlowC'] = array(
    totalC_org_high_sev*0.81)  # C
initvals['high sev burn Histosol']['uPyC'] = array(totalC_org_high_sev*0.15)
paramsets['high sev burn Histosol'] = copy.deepcopy(params)
envir_params['high sev burn Histosol'] = copy.deepcopy(envir_params)
envir_params['high sev burn Histosol']['thetamin'] = array(0.5)
envir_params['high sev burn Histosol']['thetamax'] = array(0.7)
envir_params['high sev burn Histosol']['porosity'] = array(0.4)


num_micro_pools = 0
num_micro_pools = where(SOM_init['MBC_1'] > 0,
                        num_micro_pools+1, num_micro_pools+0)
num_micro_pools = where(SOM_init['MBC_2'] > 0,
                        num_micro_pools+1, num_micro_pools+0)
num_micro_pools = where(SOM_init['MBC_3'] > 0,
                        num_micro_pools+1, num_micro_pools+0)
num_micro_pools = where(SOM_init['MBC_4'] > 0,
                        num_micro_pools+1, num_micro_pools+0)

# Should also include a function that makes sure that each simulation is using the same number of microbial pools!

# Set up a data structure to hold the results of the different simulations
results = {}
# Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
# Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
# Inputs are empty because this is running as an incubation without any constant inputs of C


### Running model with observed temperature data from Fort Smith:
# Import experimental respiration data and list of site IDs (6 of 12 sites) to be used for model parameterization
df_temp = pd.read_csv('data/daily-temp-input.csv')
# df_temp = pd.read_csv('data/static-temp-input.csv')
Tmin = df_temp['minT']
Tmax = df_temp['maxT']

# Tmin=18
# Tmax=24

### Run a simulation using the ODE solver
# for functype in initvals:
#     results[functype] = CORPSE_solvers.run_models_ODE(Tmin=Tmin,Tmax=Tmax,thetamin=envir_params[functype]['thetamin'],
#                                                       thetamax=envir_params[functype]['thetamax'],
#                                             times=t,inputs={},clay=2.5,initvals=initvals[functype],params=paramsets[functype])


### Run a simulation using the explicit iterator instead of the ODE solver
for functype in initvals:
    results[functype] = CORPSE_solvers.run_models_iterator(Tmin=Tmin,Tmax=Tmax,thetamin=envir_params[functype]['thetamin'],
                                                      thetamax=envir_params[functype]['thetamax'],
                                            times=t,inputs={},clay=2.5,initvals=initvals[functype],params=paramsets[functype])



### Code to add inputs to unburned soils:
    
# total_inputs = totalC_sandy_unburned*0.5  # 5% of total C input to system

# inputs_init = {'uFastC': total_inputs*0.4,
#                'uSlowC': total_inputs*0.6,
#                'uNecroC': total_inputs*0,
#                'uPyC': total_inputs*0}

# inputs = {}

# inputs['no burn Gleysol'] = copy.deepcopy(inputs_init)
# inputs['low sev burn Gleysol'] = copy.deepcopy(inputs_init)
# inputs['high sev burn Gleysol'] = copy.deepcopy(inputs_init)
# inputs['no burn Histosol'] = copy.deepcopy(inputs_init)
# inputs['low sev burn Histosol'] = copy.deepcopy(inputs_init)
# inputs['high sev burn Histosol'] = copy.deepcopy(inputs_init)

# inputs['low sev burn Gleysol']['uFastC'] = 0 
# inputs['low sev burn Gleysol']['uSlowC'] = 0 

# inputs['high sev burn Gleysol']['uFastC'] = 0 
# inputs['high sev burn Gleysol']['uSlowC'] = 0 

# inputs['no burn Histosol']['uFastC'] = totalC_org_unburned*0.1*0.4
# inputs['no burn Histosol']['uSlowC'] = totalC_org_unburned*0.1*0.6

# inputs['low sev burn Histosol']['uFastC'] = 0 
# inputs['low sev burn Histosol']['uSlowC'] = 0 

# inputs['high sev burn Histosol']['uFastC'] = 0 
# inputs['high sev burn Histosol']['uSlowC'] = 0 


# # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
# initvals['high sev burn Histosol'] = copy.deepcopy(SOM_init)
# # Start with low initial microbial biomass. Assumes this simulation starts right after the fire
# initvals['high sev burn Histosol']['MBC_1'] = array(
#     totalC_org_high_sev*0.0078*0.93)


# for functype in initvals:
#     results[functype] = CORPSE_solvers.run_models_iterator(Tmin=Tmin, Tmax=Tmax, thetamin=envir_params[functype]['thetamin'],
#                                                            thetamax=envir_params[functype]['thetamax'],
#                                                            times=t, inputs=dict([(k, inputs[functype][k]) for k in inputs[functype]]), clay=2.5, initvals=initvals[functype], params=paramsets[functype])
