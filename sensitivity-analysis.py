# -*- coding: utf-8 -*-
"""
Created on Fri May 10 08:33:21 2024

@author: danab
"""

## Varies each C pool and parameter related to microbial activity +/- 25% and 
#    compiles model output (cumulative C-CO2 respired over 70 days) into one
#    dataframe. This is the starting data to calculate model sensitivity (S.index)

import pandas as pd
from numpy import arange
import Whitman_sims
import copy
import CORPSE_solvers

t=arange(0,70/365,1/365) # t should be the length of the incubations that we are comparing to.
            
initvals_SA={}

paramsets_SA={}
# Create timestamp to group files later:
import time
timestr = time.strftime('%Y%m%d')


df = pd.DataFrame({"time":['time'],'simulation':['simul'],
                    'sensitivity_test':['name'],
                    'cumulative_C_CO2': ['CO2'],
                    'output_MBC':['MBC'],
                    'initial_value':['variable1'],
                    'value':[0],
                    'pool':['pool'],
                    # 'uFastC': ['uFastC'],
                    # 'uSlowC': ['uSlowC'],
                    # 'uNecroC': ['uNecroC'],
                    # 'uPyC': ['uPyC'],
                    # 'livingMicrobeC': ['livingMicrobeC']
                    })  

#### Trying to calculate S index: 
# Create list of pools to test: 
pools = ['uFastC', 'uNecroC','uSlowC', 'uPyC']

# For each simulation (burned or not burned), run separate sensitivity test
for simul in list(Whitman_sims.initvals.keys()):
     
    # For each pool (fast, slow, py, necro C), test sensitivity of model to pool size
    for i in pools:
        name = f'sensitivity_{simul}_{i}'
        initvals_SA[simul]=copy.deepcopy(Whitman_sims.initvals[simul])       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
        variable1=initvals_SA[simul][i]
        paramsets_SA[simul]=copy.deepcopy(Whitman_sims.paramsets[simul])

        for v1 in [variable1*0.75, variable1*1.0, variable1*1.25]:
            initvals_SA[simul][i]=v1
            
            
            results={}
        # Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
        # Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
        # Inputs are empty because this is running as an incubation without any constant inputs of C
    
            results[simul] = CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=0.5,thetamax=0.7,
                                                    times=t,inputs={},clay=2.5,initvals=initvals_SA[simul],params=paramsets_SA[simul])
    
    
            # Create csv file containing initial parameters and cumulative C-CO2 respired and simulation length
            # initial_pools = initvals_SA[simul]  # initial values
            # df_initial_pools = pd.DataFrame.from_dict(initial_pools, orient='index', columns = ['initial_pools_size'])
    
            # starting_params = paramsets_SA[simul]
            # df_starting_params = pd.DataFrame.from_dict(starting_params, orient='columns')
    
            output = pd.DataFrame(results[simul][0]) # Create dataframe containing simulation output
            end_time = int(t.max()*365) # Create a variable containing the last timestamp in the simulation
            output=results[simul][0].iloc[end_time]
          
            df_temp = pd.DataFrame({"time":[t],
                                    'simulation':[simul],
                                'sensitivity_test':[i],
                                'cumulative_C_CO2': [output['CO2']],
                                'output_MBC':[output['livingMicrobeC']],
                                'initial_value':[variable1],
                                'value':[v1],
                                'pool':[i],
                                # 'uFastC': ['uFastC'],
                                # 'uSlowC': ['uSlowC'],
                                # 'uNecroC': ['uNecroC'],
                                # 'uPyC': ['uPyC'],
                                # 'livingMicrobeC': ['livingMicrobeC']
                                })  
    
            df = pd.concat([df, df_temp])
            
            output = pd.DataFrame(results[simul][0]) # Create dataframe containing simulation output
            end_time = int(t.max()*365) # Create a variable containing the last timestamp in the simulation
            output=results[simul][0].iloc[end_time]
  
    
  
# Create list of parameters to test: 
parameter_pool = ['Fast','Slow','Necro','Py']
parameter_list = ['vmaxref','eup', 'kC','Ea']

# Run sensitivity tests now for all the initial variables
for simul in list(Whitman_sims.initvals.keys()):
        
    for C_pool in parameter_pool:
        for param in parameter_list:
            name = f'sensitivity_{simul}_{param}_{C_pool}'
            initvals_SA[simul]=copy.deepcopy(Whitman_sims.initvals[simul])       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation                
            paramsets_SA[simul]=copy.deepcopy(Whitman_sims.paramsets[simul])
            param_initial=paramsets_SA[simul][param][C_pool]

            for iteration in [param_initial*0.75, param_initial*1.0, param_initial*1.25]:
                paramsets_SA[simul][param][C_pool]=iteration
                
                results={}
                # Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
                # Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
                # Inputs are empty because this is running as an incubation without any constant inputs of C
        
                results[simul] = CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=0.5,thetamax=0.7,
                                                        times=t,inputs={},clay=2.5,initvals=initvals_SA[simul],params=paramsets_SA[simul])
        
        
                # Create csv file containing initial parameters and cumulative C-CO2 respired and simulation length
                # initial_pools = initvals_SA[simul]  # initial values
                # df_initial_pools = pd.DataFrame.from_dict(initial_pools, orient='index', columns = ['initial_pools_size'])
        
                # starting_params = paramsets_SA[simul]
                # df_starting_params = pd.DataFrame.from_dict(starting_params, orient='columns')
        
                output = pd.DataFrame(results[simul][0]) # Create dataframe containing simulation output
                end_time = int(t.max()*365) # Create a variable containing the last timestamp in the simulation
                output=results[simul][0].iloc[end_time]
              
                df_temp = pd.DataFrame({'simulation':[simul],
                                   'sensitivity_test':[param + "_" + C_pool],
                                   'cumulative_C_CO2': [output['CO2']],
                                   'output_MBC':[output['livingMicrobeC']],
                                   'initial_value':[param_initial],
                                   'value':[iteration],
                                   'pool':[C_pool],
                                   # 'uFastC': ['uFastC'],
                                   # 'uSlowC': ['uSlowC'],
                                   # 'uNecroC': ['uNecroC'],
                                   # 'uPyC': ['uPyC'],
                                   # 'livingMicrobeC': ['livingMicrobeC']
                                   })  
        
                df = pd.concat([df, df_temp])
                

                    
                    
df.to_csv('model_output/sensitivity_tests/' + timestr +'_SA_init_values.csv')   
   
            



            
# ### Initial sensitivity code...
#     for init_slowC in [slow_C*0.8, slow_C *0.9, slow_C, slow_C*1.1, slow_C*1.2]:
#         # for init_vmax in [slow_vmax*0.8, slow_vmax*0.9, slow_vmax, slow_vmax*1.1, slow_vmax*1.2]:
#             name=f'sensitivity_{simul}_slowC_{init_slowC:1.4f}_slowvmax_{init_vmax:1.4f}'
#             initvals_SA[simul]=copy.deepcopy(initvals[simul])       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
#             initvals_SA[simul]['uSlowC']=init_slowC
#             paramsets_SA[simul]=copy.deepcopy(paramsets[simul])
#             paramsets_SA[simul]['vmaxref']['Slow']=init_vmax    
#             #paramsets_SA['sensitivity']['eup']['Slow'] = sloweup  
