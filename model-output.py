# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:36:32 2024

@author: danab
"""
## Generate data file with model output from all CORPSE simulations to use in R

import pandas as pd
from numpy import arange
import Whitman_sims
import CORPSE_solvers
import CORPSE_array

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
                   'livingMicrobeC': [0]
                    })  

          

# For each simulation (burned or not burned), run separate sensitivity test
for simul in list(Whitman_sims.initvals.keys()):
    results={}
    # Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
    # Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
    # Inputs are empty because this is running as an incubation without any constant inputs of C
    
    results[simul] = CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=0.5,thetamax=0.7,
                                                    times=t,inputs={},clay=2.5,initvals=Whitman_sims.initvals[simul],params=Whitman_sims.paramsets[simul])
    
    # Ned to include cumulative C-CO2 as percent of initial total C
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[simul][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[simul][0], 'p')
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
            'livingMicrobeC': output['livingMicrobeC']
              }
    
    df_temp= {}
    df_temp = pd.DataFrame(temp)
    df_temp.loc[:, "time"]=t
    df_temp.loc[:, "simulation"]=simul
    df_temp.loc[:, "total_initial_C"]=total_initial_C
    
      
    # df = pd.concat([df, df_temp])







    df = pd.concat([df, df_temp])

    
    
    

# df_initial_pools.to_csv('model_output/sensitivity_tests/' + timestr + '_SA_initial_pools.csv')
# df_starting_params.to_csv('model_output/sensitivity_tests/' + timestr + '_SA_initial_params.csv')
df.to_csv('model_output/' + timestr + '_model-output.csv')








