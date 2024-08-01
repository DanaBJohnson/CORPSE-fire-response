# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:36:32 2024

@author: danab
"""
## Generate data file with model output from all CORPSE simulations to use in R

import pandas as pd
from numpy import arange
import Microbial_sims
import Microbial_CORPSE_solvers
import Microbial_CORPSE_array

t=arange(0,70/365,1/365) # t should be the length of the incubations that we are comparing to.
            
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
    
    results[simul] = Microbial_CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=0.5,thetamax=0.7,
                                                    times=t,inputs={},clay=2.5,initvals=Microbial_sims.initvals[simul],params=Microbial_sims.paramsets[simul])
    
    # Ned to include cumulative C-CO2 as percent of initial total C
    #totalC=Microbial_CORPSE_array.sumCtypes(Microbial_sims.results[simul][0], 'u')+Microbial_CORPSE_array.sumCtypes(Microbial_sims.results[simul][0], 'p')
    #total_initial_C = totalC[0]
    totalC=Microbial_CORPSE_solvers.totalCarbon(Microbial_sims.results[simul][0], Microbial_CORPSE_array.microbial_pools)
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
df.to_csv('microbial_model_output/' + timestr + '_model-output.csv')








