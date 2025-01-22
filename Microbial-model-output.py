# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:36:32 2024

@author: danab
"""
## Generate data file with model output from all CORPSE simulations to use in R

import pandas as pd
from numpy import arange

import Microbial_sims_single_MBC_pool as Microbial_sims
# If using Microbial_sims_single_MBC_pool
pools="one"

# import Microbial_sims as Microbial_sims
# # If using Microbial_sims
# pools="two"


import CORPSE_solvers
import CORPSE_array

t=arange(0,365/365,1/365) # t should be the length of the incubations that we are comparing to.
            
initvals_SA={}

paramsets_SA={}
# Create timestamp to group files later:
import time
timestr = time.strftime('%Y%m%d')


df = pd.DataFrame({"time":[0],
                   'simulation':['simul'],
                   'cumulative_C_CO2': [0],
                   'total_initial_C':[0],
                   'uFastC': [0],
                   'uSlowC': [0],
                   'uNecroC': [0],
                   'uPyC': [0],
                   'MBC_1': [0],
                   'MBC_2':[0]
                    })  



# For each simulation (burned or not burned), run separate sensitivity test
for simul in list(Microbial_sims.initvals.keys()):
    results={}
    # Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
    # Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
    # Inputs are empty because this is running as an incubation without any constant inputs of C
    
    results[simul] = CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=0.5,thetamax=0.7,
                                                    times=t,inputs={},clay=2.5,initvals=Microbial_sims.initvals[simul],params=Microbial_sims.paramsets[simul])
    
    # Ned to include cumulative C-CO2 as percent of initial total C
    #totalC=CORPSE_array.sumCtypes(Microbial_sims.results[simul][0], 'u')+CORPSE_array.sumCtypes(Microbial_sims.results[simul][0], 'p')
    #total_initial_C = totalC[0]
    totalC=CORPSE_solvers.totalCarbon(Microbial_sims.results[simul][0], CORPSE_array.microbial_pools)
    total_initial_C = totalC[0]
    
    
    output = pd.DataFrame(results[simul][0]) # Create dataframe containing simulation output
    end_time = int(t.max()*365) # Create a variable containing the last timestamp in the simulation
  
    # df_temp = pd.DataFrame({"time":[t],
    #                         'simulation':[simul],
    #                         'cumulative_C_CO2': [output['CO2']],
    #                         'livingMicrobeC':[output['livingMicrobeC']],
    #                         'uFastC': [output['uFastC']],
    #                         'uSlowC': [output['uSlowC']],
    #                         'uNecroC': [output['uNecroC']],
    #                         'uPyC': [output['uPyC']],
    #                         })  

    temp = {}
    temp = {'cumulative_C_CO2': output['CO2'],
            'uFastC': output['uFastC'],
            'uSlowC': output['uSlowC'],
            'uNecroC': output['uNecroC'],
            'uPyC': output['uPyC'],
            'MBC_1': output['MBC_1'],
            'MBC_2': output['MBC_2']
              }
    
    df_temp= {}
    df_temp = pd.DataFrame(temp)
    df_temp.loc[:, "time"]=t
    df_temp.loc[:, "simulation"]=simul
    df_temp.loc[:, "total_initial_C"]=total_initial_C

    df = pd.concat([df, df_temp])

   
# df_initial_pools.to_csv('model_output/sensitivity_tests/' + timestr + '_SA_initial_pools.csv')
# df_starting_params.to_csv('model_output/sensitivity_tests/' + timestr + '_SA_initial_params.csv')
df.to_csv('microbial_model_output/' + timestr + '-' + pools + '_pool_model-output.csv')






# Rather ugly code to organize and safe initial parameters for each run:
df_vmax = pd.DataFrame.from_dict(data=Microbial_sims.params['vmaxref'], orient='index') 
df_vmax['param']='vmaxref'
df_vmax['MBC'] = ['MBC_1','MBC_2','MBC_3','MBC_4']

df_kC = pd.DataFrame.from_dict(data=Microbial_sims.params['kC'], orient='index') 
df_kC['param']='kC'
df_kC['MBC'] = ['MBC_1','MBC_2','MBC_3','MBC_4']

df_eup = pd.DataFrame.from_dict(data=Microbial_sims.params['eup'], orient='index') 
df_eup['param']='eup'
df_eup['MBC'] = ['MBC_1','MBC_2','MBC_3','MBC_4']

df_params1 = pd.concat([df_vmax,df_kC, df_eup])
df_params1['No pool parameter value'] = 'NA'
df_params1 = df_params1.reset_index(drop=True)

df_params2 = pd.DataFrame([['Ea', Microbial_sims.params['Ea']['Fast'], Microbial_sims.params['Ea']['Slow'], 
                    Microbial_sims.params['Ea']['Py'], Microbial_sims.params['Ea']['Necro']],
                   ['protection_rate',Microbial_sims.params['protection_rate']['Fast'], Microbial_sims.params['protection_rate']['Slow'],
                    Microbial_sims.params['protection_rate']['Py'], Microbial_sims.params['protection_rate']['Necro']]],
                  columns = ['param','Fast', 'Slow','Py','Necro'])
df_params2['MBC'] = 'NA'
df_params2['No pool parameter value'] = 'NA'


df_Tmic = pd.DataFrame.from_dict(data=Microbial_sims.params['Tmic'], orient='index')
df_Tmic['param'] = 'Tmic'
df_Tmic['MBC'] = ['MBC_1','MBC_2','MBC_3','MBC_4']

df_minMic = pd.DataFrame.from_dict(data=Microbial_sims.params['minMicrobeC'], orient='index')
df_minMic['param'] = 'minMicrobeC'
df_minMic['MBC'] = ['MBC_1','MBC_2','MBC_3','MBC_4']

df_et = pd.DataFrame.from_dict(data=Microbial_sims.params['et'], orient='index')
df_et['param'] = 'et'
df_et['MBC'] = ['MBC_1','MBC_2','MBC_3','MBC_4']


data = [['gas_diffusion_exp', 0.6],
        ['substrate_diffusion_exp', 1.5],
        ['tProtected', 75.0], 
        ['new_resp_units', True]]

df_others = pd.DataFrame(data, columns = ['param',0])

df_params3 = pd.concat([df_Tmic, df_minMic, df_et, df_others])
df_params3 = df_params3.rename(columns={0: "No pool parameter value"})
df_params3 = df_params3.reset_index(drop=True)


df_params3['Fast'] = 'NA'
df_params3['Slow'] = 'NA'
df_params3['Py'] = 'NA'
df_params3['Necro'] = 'NA'

df_params = pd.concat([df_params1, df_params2, df_params3])
df_params['Model version'] = pools

df_params.to_csv('microbial_model_output/' + timestr + '-' + pools + '_pool_parameters.csv')
