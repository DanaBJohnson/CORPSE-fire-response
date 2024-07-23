# -*- coding: utf-8 -*-
"""
Created on Fri May 10 13:53:28 2024

@author: danab
"""
## Varies C pools and parameters related to microbial activity +/- 10 and 20% and 
#    compiles the full model output (C in each pool varying with time over the 
#    course of the simulation) into one dataframe. This dataframe is saved to 
#    the "calibrations" folder and used for calibration and parameterization. 

# Focus on varying parameters that the model is most sensitive to:
    # Unburned:
        # uFastC
        # uSlowC
        # uNecroC
        # vmaxref_Slow
        # kC_Slow
    # Burned:
        # uFastC
        # uSlowC
        # vmaxref_Fast
        # vmaxref_Slow
        # eup_fast



import pandas as pd
from numpy import arange
import Whitman_sims
import copy
import CORPSE_solvers

# Setting up sensitivity analyses: 
df = pd.DataFrame({'time':['time'],'simulation':['simul'],
                    'calibration_test':['name'],
                    'cumulative_C_CO2': ['CO2'],
                    # 'uFastC': ['uFastC'],
                    # 'uSlowC': ['uSlowC'],
                    # 'uNecroC': ['uNecroC'],
                    # 'uPyC': ['uPyC'],
                    # 'livingMicrobeC': ['livingMicrobeC']
                    'Value1_initial':['Value1'],
                    'Value2_initial':['Value2'],
                    'Value3_initial':['Value3'],
                    'Parameter1_initial':['Param1'],
                    'Parameter2_initial':['Param2'],
                    'Parameter3_initial':['Param3']
                    })  
  
t=arange(0,70/365,1/365) # t should be the length of the incubations that we are comparing to.
            
initvals_SA={}

paramsets_SA={}
# Create timestamp to group files later:
import time
timestr = time.strftime('%Y%m%d')

# Effect of varying CUE?
for simul in list(Whitman_sims.initvals.keys()):
    
    # Combinations to test:
    # combos = {"pool": ['uFastC', 'uFastC', 'uSlowC','uSlowC', 'uPyC', 'uPyC'],
    #           "parameter_pool": ['Fast', 'Fast', 'Slow','Slow', 'Py', 'Py'],
    #           "parameter": ['vmaxref','eup', 'vmaxref','eup', 'vmaxref','eup']}

    # combos = {"pool": ['uFastC', 'uFastC', 'uSlowC','uSlowC'],
    #           "parameter_pool": ['Fast', 'Fast', 'Slow','Slow'],
    #           "parameter": ['vmaxref','eup', 'vmaxref','eup']}
       
    # combos = {"pool": ['uFastC', 'uFastC'],
    #           "parameter_pool": ['Fast', 'Slow'],
    #           "parameter": ['vmaxref','vmaxref']}
    
    
    # simul = 'no burn sandy soil'
    # i = 0
    # v1=variable1*0.9
    # v2=variable2*0.9
    
    # for i in list(range(0,len(combos['pool']))):
    #     name1=combos['pool'][i]
    #     name2=(combos['pool'][i] + '' + combos['parameter'][i])
    #     initvals_SA[simul]=copy.deepcopy(Whitman_sims.initvals[simul])       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
    #     paramsets_SA[simul]=copy.deepcopy(Whitman_sims.paramsets[simul])
    #     variable1=initvals_SA[simul][name1]
    #     variable2=paramsets_SA[simul][combos['parameter'][i]][combos['parameter_pool'][i]]

       
    name1='uFastC'
    name2='uSlowC'
    name3='uNecroC'
    name4='vmaxref'
    pool4='Fast'
    name5='vmaxref'
    pool5='Slow'
    name6='eup'
    pool6='Fast'
    
    initvals_SA[simul]=copy.deepcopy(Whitman_sims.initvals[simul])       # Makes a copy of the default initial values. Need to use deepcopy so changing the value here doesn't change it for every simulation
    paramsets_SA[simul]=copy.deepcopy(Whitman_sims.paramsets[simul])
    variable1=initvals_SA[simul][name1]
    variable2=initvals_SA[simul][name2]
    variable3=initvals_SA[simul][name3]
    
    # add eup_Fast, Vmax_F, and Vmax_S
    parameter1=paramsets_SA[simul][name4][pool4]
    parameter2=paramsets_SA[simul][name5][pool5]
    parameter3=paramsets_SA[simul][name6][pool6]
    
    # for v1 in [variable1*0.8,variable1*0.9, variable1*1.0, variable1*1.1,variable1*1.2]:
        # for v2 in [variable2*0.8,variable2*0.9, variable2*1.0, variable2*1.1, variable2*1.2]:
            # for v3 in [variable3*0.8,variable3*0.9, variable3*1.0, variable3*1.1,variable3*1.2]:
    for p1 in [parameter1*0.8,parameter1*0.9, parameter1, parameter1*1.1,parameter1*1.2]:
        for p2 in [parameter2*0.8,parameter2*0.9, parameter2, parameter2*1.1,parameter2*1.2]:
            for p3 in [parameter3*0.8,parameter3*0.9, parameter3, parameter3*1.1,parameter3*1.2]:
                
                
                            # name=f'calibration_{simul}_{name1}_{v1:1.4f}_{name2}_{v2:1.4f}_{name3}_{v3:1.4f}_{name4}_{pool4}_{p1:1.4f}_{name5}_{pool5}_{p2:1.4f}_{name6}_{pool6}_{p3:1.4f}'                           
                            # name=f'calibration_{simul}_{name1}_{v1:1.4f}_{name2}_{v2:1.4f}_{name3}_{v3:1.4f}'
                            # initvals_SA[simul][name1]=v1
                            # initvals_SA[simul][name2]=v2
                            # initvals_SA[simul][name3]=v3
                            name=f'calibration_{simul}_{name4}{pool4}_{p1:1.4f}_{name5}{pool5}_{p2:1.4f}_{name6}{pool6}_{p3:1.4f}'
                            paramsets_SA[simul][name4][pool4]=p1
                            paramsets_SA[simul][name5][pool5]=p2
                            paramsets_SA[simul][name6][pool6]=p3
                            
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
                                # end_time = int(t.max()*365) # Create a variable containing the last timestamp in the simulation
                                # output=results[simul][0].iloc[end_time]
                            temp = {}
                            temp = {'cumulative_C_CO2': output['CO2'],
                                        # 'uFastC': output['uFastC'],
                                        # 'uSlowC': output['uSlowC'],
                                      # 'uNecroC': output['uNecroC'],
                                      # 'uPyC': output['uPyC'],
                                      # 'livingMicrobeC': output['livingMicrobeC']
                                      }
                            df_temp= {}
                            df_temp = pd.DataFrame(temp)
                            df_temp.loc[:, "simulation"]=simul
                            df_temp.loc[:, "calibration_test"]=name
                            df_temp.loc[:, "time"]=t
                            df_temp.loc[:, "Value1_initial"]=variable1
                            df_temp.loc[:, "Value2_initial"]=variable2
                            df_temp.loc[:, "Value3_initial"]=variable3
                            df_temp.loc[:, 'Parameter1_initial']=parameter1
                            df_temp.loc[:, 'Parameter2_initial']=parameter2
                            df_temp.loc[:, 'Parameter3_initial']=parameter3
                            df = pd.concat([df, df_temp])
                            


        # for v1 in [variable1*0.5, variable1*0.6, variable1*0.7, variable1*0.8, variable1*0.9, variable1*1.0, variable1*1.1, variable1*1.2, variable1*1.3,variable1*1.4, variable1*1.5]:
        #     for v2 in [variable2*0.8, variable2*0.9, variable2*1.0, variable2*1.1, variable2*1.2]:
        #         name=f'calibration_{simul}_{name1}_{v1:1.4f}_{name2}_{v2:1.4f}_{name3}_{v3:1.4f}'
        #         name=f'calibration_{simul}_{name1}_{v1:1.4f}'
        #         initvals_SA[simul][name1]=v1
        #             #paramsets_SA[simul][combos['parameter'][i]][combos['parameter_pool'][i]]=v2  
                
              
        #         # Set up a data structure to hold the results of the different simulations
        #         results={}
        #         # Goes through each functional type and runs a simulation using the appropriate set of parameters and initial values
        #         # Simulations are assuming a constant temperature of 20 C and constant moisture of 60% of saturation
        #         # Inputs are empty because this is running as an incubation without any constant inputs of C
        
        #         results[simul] = CORPSE_solvers.run_models_ODE(Tmin=18.0,Tmax=24.0,thetamin=0.5,thetamax=0.7,
        #                                                     times=t,inputs={},clay=2.5,initvals=initvals_SA[simul],params=paramsets_SA[simul])
        
        
        #             # Create csv file containing initial parameters and cumulative C-CO2 respired and simulation length
        #             # initial_pools = initvals_SA[simul]  # initial values
        #             # df_initial_pools = pd.DataFrame.from_dict(initial_pools, orient='index', columns = ['initial_pools_size'])
        
        #             # starting_params = paramsets_SA[simul]
        #             # df_starting_params = pd.DataFrame.from_dict(starting_params, orient='columns')
        
        #         output = pd.DataFrame(results[simul][0]) # Create dataframe containing simulation output
        #             # end_time = int(t.max()*365) # Create a variable containing the last timestamp in the simulation
        #             # output=results[simul][0].iloc[end_time]
        #         temp = {}
        #         temp = {'cumulative_C_CO2': output['CO2'],
        #                     # 'uFastC': output['uFastC'],
        #                     # 'uSlowC': output['uSlowC'],
        #                   # 'uNecroC': output['uNecroC'],
        #                   # 'uPyC': output['uPyC'],
        #                   # 'livingMicrobeC': output['livingMicrobeC']
        #                   }
        #         df_temp= {}
        #         df_temp = pd.DataFrame(temp)
        #         df_temp.loc[:, "simulation"]=simul
        #         df_temp.loc[:, "calibration_test"]=name
        #         df_temp.loc[:, "time"]=t
        
        #         df = pd.concat([df, df_temp])
                

# df_initial_pools.to_csv('model_output/sensitivity_tests/' + timestr + '_Calibration_initial_pools.csv')
# df_starting_params.to_csv('model_output/sensitivity_tests/' + timestr + '_Calibration_initial_params.csv')
# df.to_csv('model_output/sensitivity_tests/' + timestr +'_Calibration_results_' + 'slowC_vs_fvmax' +'.csv')
df.to_csv('model_output/calibrations/' + timestr +'_Calibration_results.csv')   


