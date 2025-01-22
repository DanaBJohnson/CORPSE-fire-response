# Import respiration data
import pandas as pd

# Import experimental respiration data and list of site IDs (6 of 12 sites) to be used for model parameterization
df = pd.read_csv('../../data/incubations/multiplexer/processed-respiration-data.csv')
sites = pd.read_csv('data/site-ID-to-use-for-calibrations.csv')
IDs = sites['site.id']
df=df[df["site.id"].isin(IDs)]


from matplotlib import pyplot
import CORPSE_array as CORPSE_array

# Import CORPSE simulation results.

# This section plots the results
# Each set of results should have the same set of pools as the initial values structure from the beginning of the simulation

################### 
# Indicate which simulation to import: either Microbial_sims (multi-MBC pool simulation) or Microbial_sims_single_MBC_pool (single-MBC pool simulation)
import Microbial_sims as Microbial_sims
# import Microbial_sims_single_MBC_pool as Microbial_sims


nrows=int(Microbial_sims.num_micro_pools)

# if (Microbial_sims.num_micro_pools)>1:
#     #nrows = len(CORPSE_array.microbial_pools)
#     fig,ax=pyplot.subplots(nrows=nrows,ncols=1,clear=True,num='CORPSE results')
#     p=0
#     c=0
#     for m in CORPSE_array.microbial_pools[0:nrows]:
#         for sim in Microbial_sims.results:
#             sim_color = ['palegreen','khaki','salmon','forestgreen','darkgoldenrod','darkred',
#                          'palegreen','khaki','salmon','forestgreen','darkgoldenrod','darkred']
#             totalC=CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'p')
#             ax[p].plot(Microbial_sims.t*365,Microbial_sims.results[sim][0][m]/totalC[0]*100, label = sim, color = sim_color[c], linewidth=3) # Cumulative % of initial C respired)
#             c=c+1
#         ax[p].set_title(m, y=1, pad=-15)
#         ax[p].set_ylabel('Microbial pool biomass \n(% initial C)')
#         p=p+1
#     ax[p-1].set_xlabel('Time (days)')
#     ax[p-1].legend(fontsize='small')
# else:
#     fig,ax=pyplot.subplots(nrows=1,ncols=1,clear=True,num='CORPSE results')
#     m='MBC_1'
#     for sim in Microbial_sims.results:
#         totalC=CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'p')
#         ax.plot(Microbial_sims.t*365,Microbial_sims.results[sim][0][m]/totalC[0]*100, label = sim) # Cumulative % of initial C respired)
#         ax.set_title(m, y=1, pad=-15)
#     ax.set_ylabel('Microbial pool biomass \n(% initial C)')
#     ax.set_xlabel('Time (days)')
#     ax.legend(fontsize='small')

# pyplot.subplots_adjust(hspace = 0.2)

# pyplot.show(block=False)
#pyplot.draw()





# # Create figure XXX
# fig,ax=pyplot.subplots(nrows=2,ncols=2,clear=True,num='CORPSE results')
# for sim in Microbial_sims.results:
#     totalC=CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'p')
#     ax[0,0].plot(Microbial_sims.t*365, Microbial_sims.results[sim][0]['uFastC']/totalC[0]*100, label = sim)
#     #ax.plot(Microbial_sims.t*365, Microbial_sims.results[sim][0]['uSlowC']/totalC[0]*100, color = 'blue', label = 'uSlowC')
#     ax[1,0].plot(Microbial_sims.t*365, Microbial_sims.results[sim][0]['uNecroC']/totalC[0]*100, label = sim)
#     ax[0,1].plot(Microbial_sims.t*365, Microbial_sims.results[sim][0]['uPyC']/totalC[0]*100, label = sim)
#     ax[1,1].plot(Microbial_sims.t*365, Microbial_sims.results[sim][0]['uSlowC']/totalC[0]*100, label = sim)

# #ax[0,0].set_xlabel('Time (days)')
# ax[0,0].set_ylabel('Percent of total C')
# ax[0,0].legend(loc='center left', prop={'size':6})
# ax[0,0].set_title('uFastC', y=1, pad=10)
# ax[1,0].set_xlabel('Time (days)')
# ax[1,0].set_ylabel('Percent of total C')
# #ax[1,0].legend(fontsize='small')
# ax[1,0].set_title('uNecroC', y=1, pad=-10)
# ax[0,1].set_xlabel('Time (days)')
# #ax[0,1].set_ylabel('Percent of total C')
# #ax[0,1].legend(fontsize='small')
# ax[0,1].set_title('uPyC')
# ax[1,1].set_xlabel('Time (days)')
# #ax[1,1].set_ylabel('Percent of total C')
# #ax[1,1].legend(fontsize='small')
# ax[1,1].set_title('uSlowC', y=1, pad=-10)

# pyplot.subplots_adjust(hspace = 0.2)

# pyplot.show(block=False)
# #pyplot.draw()



# Put the model and experimental results on the same figure.
fig,ax=pyplot.subplots(nrows=2,ncols=1,clear=True)
for sim in Microbial_sims.results:
    totalC=CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'p')
    if (sim=='low sev burn Gleysol' or sim=='high sev burn Gleysol'): 
        ax[0].plot(Microbial_sims.t*365,Microbial_sims.results[sim][0]['CO2']/totalC[0]*100, label = 'Model results: '+sim, linewidth=3.0) # Cumulative % of initial C respired
    elif (sim=='low sev burn Histosol' or sim=='high sev burn Histosol'): 
        ax[1].plot(Microbial_sims.t*365,Microbial_sims.results[sim][0]['CO2']/totalC[0]*100, label = 'Model results: '+sim, linewidth=3.0) # Cumulative % of initial C respired


# Add laboratory respiration data to figure:
# no_burn_Pinus = df[(df["burn.trtmt"]=='control') & (df['vegetation']=='Pinus_banksiana')]
# x_noburn_Pinus = no_burn_Pinus['whole.days.since.wet.up']
# y_noburn_Pinus = no_burn_Pinus['new.cum_CO2C_g']/no_burn_Pinus['initial.total.C.g']*100
burn_30s_Pinus = df[(df["burn.trtmt"]=='low severity') & (df['vegetation']=='Pinus_banksiana')]
x_30burn_Pinus = burn_30s_Pinus['whole.days.since.wet.up']
y_30burn_Pinus = burn_30s_Pinus['new.cum_CO2C_g']/burn_30s_Pinus['initial.total.C.g']*100
burn_120s_Pinus = df[(df["burn.trtmt"]=='high severity') & (df['vegetation']=='Pinus_banksiana')]
x_120burn_Pinus = burn_120s_Pinus['whole.days.since.wet.up']
y_120burn_Pinus = burn_120s_Pinus['new.cum_CO2C_g']/burn_120s_Pinus['initial.total.C.g']*100
# no_burn_Picea = df[(df["burn.trtmt"]=='control') & (df['vegetation']=='Picea_spp.')]
# x_noburn_Picea = no_burn_Picea['whole.days.since.wet.up']
# y_noburn_Picea = no_burn_Picea['new.cum_CO2C_g']/no_burn_Picea['initial.total.C.g']*100
burn_30s_Picea = df[(df["burn.trtmt"]=='low severity') & (df['vegetation']=='Picea_spp.')]
x_30burn_Picea = burn_30s_Picea['whole.days.since.wet.up']
y_30burn_Picea = burn_30s_Picea['new.cum_CO2C_g']/burn_30s_Picea['initial.total.C.g']*100
burn_120s_Picea = df[(df["burn.trtmt"]=='high severity') & (df['vegetation']=='Picea_spp.')]
x_120burn_Picea = burn_120s_Picea['whole.days.since.wet.up']
y_120burn_Picea = burn_120s_Picea['new.cum_CO2C_g']/burn_120s_Picea['initial.total.C.g']*100

# ax[0].scatter(x_noburn_Pinus, y_noburn_Pinus, s=2, color = 'cornflowerblue', label = "Lab results: No burn")
ax[0].scatter(x_30burn_Pinus, y_30burn_Pinus, s=2, color = 'cornflowerblue', label = "Lab results: low severity burn")
ax[0].scatter(x_120burn_Pinus, y_120burn_Pinus, s=2, color = 'sandybrown', label = "Lab results: high sev burn ")
# ax[1].scatter(x_noburn_Picea, y_noburn_Picea, s=2, color = 'cornflowerblue', label = "Lab results: No burn")
ax[1].scatter(x_30burn_Picea, y_30burn_Picea, s=2, color = 'cornflowerblue', label = "Lab results: low severity burn")
ax[1].scatter(x_120burn_Picea, y_120burn_Picea, s=2, color = 'sandybrown', label = "Lab results: high sev burn")

# Format axes and legend
ax[0].set_xlabel('Time (days)')
ax[0].set_ylabel('Cumulative C-CO$_2$ respired\n(% of initial C)')
ax[0].set_title('Multi-pool model: Gleysol soil', y=1, pad=-15)
# ax[0].set_xlim([0,35])
#ax[0].legend(loc = 'upper left', prop={'size':8})
ax[1].set_xlabel('Time (days)')
ax[1].set_ylabel('Cumulative C-CO$_2$ respired\n(% of initial C)')
ax[1].set_title('Histosol soil', y=1, pad=-15)
ax[1].legend(loc = 'upper left', prop={'size':8})
# ax[1].set_xlim([0,35])


pyplot.show() # Display plot
   
    
   
# fig,ax=pyplot.subplots(nrows=1,ncols=1,clear=True,num='Vmax')
# ax.scatter(t*365, Vmax['Fast'], color = 'green', s=3)
# ax.scatter(t*365, Vmax['Slow'], color = 'red',s=3)
# ax.scatter(t*365, Vmax['Necro'], color = 'brown',s=3)
# ax.scatter(t*365, Vmax['Py'], color='black',s=3)
# # ax.scatter(t*365, Tmin, color='blue',s=2)
# # ax.scatter(t*365, Tmax, color='red',s=2)
# # ax.scatter(t*365, T, color='green',s=5)
# pyplot.show() # Display plot


   
    
# Put the model and experimental results on the same figure.
fig,ax=pyplot.subplots(nrows=2,ncols=1,clear=True)
for sim in Microbial_sims.results:
    totalC=CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'p')
    if (sim=='no burn Gleysol'): 
        ax[0].plot(Microbial_sims.t*365,Microbial_sims.results[sim][0]['CO2']/totalC[0]*100, label = 'Model results: '+sim, linewidth=3.0) # Cumulative % of initial C respired
    elif (sim=='no burn Histosol'):
        ax[1].plot(Microbial_sims.t*365,Microbial_sims.results[sim][0]['CO2']/totalC[0]*100, label = 'Model results: '+sim, linewidth=3.0) # Cumulative % of initial C respired


# Add laboratory respiration data to figure:
no_burn_Pinus = df[(df["burn.trtmt"]=='control') & (df['vegetation']=='Pinus_banksiana')]
x_noburn_Pinus = no_burn_Pinus['whole.days.since.wet.up']
y_noburn_Pinus = no_burn_Pinus['new.cum_CO2C_g']/no_burn_Pinus['initial.total.C.g']*100
# burn_30s_Pinus = df[(df["burn.trtmt"]=='low severity') & (df['vegetation']=='Pinus_banksiana')]
# x_30burn_Pinus = burn_30s_Pinus['whole.days.since.wet.up']
# y_30burn_Pinus = burn_30s_Pinus['new.cum_CO2C_g']/burn_30s_Pinus['initial.total.C.g']*100
# burn_120s_Pinus = df[(df["burn.trtmt"]=='high severity') & (df['vegetation']=='Pinus_banksiana')]
# x_120burn_Pinus = burn_120s_Pinus['whole.days.since.wet.up']
# y_120burn_Pinus = burn_120s_Pinus['new.cum_CO2C_g']/burn_120s_Pinus['initial.total.C.g']*100
no_burn_Picea = df[(df["burn.trtmt"]=='control') & (df['vegetation']=='Picea_spp.')]
x_noburn_Picea = no_burn_Picea['whole.days.since.wet.up']
y_noburn_Picea = no_burn_Picea['new.cum_CO2C_g']/no_burn_Picea['initial.total.C.g']*100
# burn_30s_Picea = df[(df["burn.trtmt"]=='low severity') & (df['vegetation']=='Picea_spp.')]
# x_30burn_Picea = burn_30s_Picea['whole.days.since.wet.up']
# y_30burn_Picea = burn_30s_Picea['new.cum_CO2C_g']/burn_30s_Picea['initial.total.C.g']*100
# burn_120s_Picea = df[(df["burn.trtmt"]=='high severity') & (df['vegetation']=='Picea_spp.')]
# x_120burn_Picea = burn_120s_Picea['whole.days.since.wet.up']
# y_120burn_Picea = burn_120s_Picea['new.cum_CO2C_g']/burn_120s_Picea['initial.total.C.g']*100

ax[0].scatter(x_noburn_Pinus, y_noburn_Pinus, s=2, color = 'cornflowerblue', label = "Lab results: No burn")
# ax[0].scatter(x_30burn_Pinus, y_30burn_Pinus, s=2, color = 'sandybrown', label = "Lab results: low severity burn")
# ax[0].scatter(x_120burn_Pinus, y_120burn_Pinus, s=2, color = 'darkseagreen', label = "Lab results: high sev burn ")
ax[1].scatter(x_noburn_Picea, y_noburn_Picea, s=2, color = 'cornflowerblue', label = "Lab results: No burn")
# ax[1].scatter(x_30burn_Picea, y_30burn_Picea, s=2, color = 'sandybrown', label = "Lab results: low severity burn")
# ax[1].scatter(x_120burn_Picea, y_120burn_Picea, s=2, color = 'darkseagreen', label = "Lab results: high sev burn")

# Format axes and legend
ax[0].set_xlabel('Time (days)')
ax[0].set_ylabel('Cumulative C-CO$_2$ respired\n(% of initial C)')
ax[0].set_title('No burn: Gleysol soil', y=1, pad=-15)
# ax[0].set_xlim([0,35])
#ax[0].legend(loc = 'upper left', prop={'size':8})
ax[1].set_xlabel('Time (days)')
ax[1].set_ylabel('Cumulative C-CO$_2$ respired\n(% of initial C)')
ax[1].set_title('No burn: Histosol soil', y=1, pad=-15)
ax[1].legend(loc = 'lower right', prop={'size':8})
# ax[1].set_xlim([0,35])


pyplot.show() # Display plot




# Check model output by summing totalC and CO2
#fig,ax=pyplot.subplots(nrows=1)

# for sim in Whitman_sims.results:
#     totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
#     ax.plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['CO2'] + totalC + Whitman_sims.results[sim][0]['MBC'])

# ax.set_xlabel('Time (days)')
# ax.set_ylabel('Total C in model')
#ax.set_title('Laboratory incubation results', y=1, pad=-15)

#pyplot.show()


# Save results file
#df_output = Whitman_sims.results['burn'][0]
#df_output.to_csv('model_output/No_burn_Picea.csv', index_label='Time')

