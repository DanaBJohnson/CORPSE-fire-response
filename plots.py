# Import respiration data
import pandas as pd

# Import experimental respiration data and list of site IDs (6 of 12 sites) to be used for model parameterization
df = pd.read_csv('../../data/incubations/multiplexer/processed-respiration-data.csv')
sites = pd.read_csv('data/site-ID-to-use-for-calibrations.csv')
IDs = sites['site.id']
df=df[df["site.id"].isin(IDs)]

    

# Plot experimental respiration data from Jack pine sites following the 'No burn' treatment. 
from matplotlib import pyplot

#df=df[df['vegetation']=='Picea_spp.'] 

fig,ax=pyplot.subplots(nrows=1)
no_burn = df[df["burn.trtmt"]=='control']
x_noburn = no_burn['whole.days.since.wet.up']
y_noburn = no_burn['g.CO2C.per.initial.g.C.per.day']*100

burn_30s = df[df["burn.trtmt"]=='low severity']
x_30burn = burn_30s['whole.days.since.wet.up']
y_30burn = burn_30s['g.CO2C.per.initial.g.C.per.day']*100

burn_120s = df[df["burn.trtmt"]=='high severity']
x_120burn = burn_120s['whole.days.since.wet.up']
y_120burn = burn_120s['g.CO2C.per.initial.g.C.per.day']*100

ax.scatter(x_noburn, y_noburn, s=3, color = 'grey', label = "no burn")
ax.scatter(x_30burn, y_30burn, s=3, color = 'orange', label = "low severity burn")
ax.scatter(x_120burn, y_120burn, s=3, color = 'red', label = "high severity burn")

ax.set_xlabel('Time (days)')
ax.set_ylabel('CO$_2$ flux rate \n(% initial C/day)')
ax.set_title('Laboratory incubation results', y=1, pad=-15)
ax.legend(loc = 'center right')

pyplot.show()


# Import CORPSE simulation results.

# This section plots the results
# Each set of results should have the same set of pools as the initial values structure from the beginning of the simulation
import Whitman_sims
import CORPSE_array

fig,ax=pyplot.subplots(nrows=2,ncols=1,clear=True,num='CORPSE results')

for sim in Whitman_sims.results:
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    #ax[0].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2'].diff()/totalC[0]*100,label=sim) # % of initial C respired each day
    ax[0].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2']/totalC[0]*100, label = sim) # Cumulative % of initial C respired
    # ax[1].plot(t*365,results[sim][0]['uFastC'],label='Simple')
    # ax[1].plot(t*365,results[sim][0]['uSlowC'],label='Complex')
    # ax[1].plot(t*365,results[sim][0]['uNecroC'],label='Necromass')
    ax[1].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['livingMicrobeC']/totalC[0]*100)
    #ax[2].plot(t*365,totalC) # How to plot total C in all four functional groups combined?

ax[0].set_xlabel('Time (days)')
ax[1].set_xlabel('Time (days)')
# ax[2].set_xlabel('Time (days)')
ax[0].set_ylabel('CO$_2$ flux rate \n(% initial C/day)')
ax[1].set_ylabel('Microbial biomass \n(% initial C)')
# ax[1].set_ylabel('SOM pools')
# ax[1].legend()
ax[0].legend(fontsize='small')
ax[0].set_title('CO$_2$ fluxes', y=1, pad=-15)
ax[1].set_title('Microbial biomass', y=1, pad=-15)

pyplot.subplots_adjust(hspace = 0.2)

pyplot.show()


# Create figure XXX
fig,ax=pyplot.subplots(nrows=2,ncols=2,clear=True,num='CORPSE results')
for sim in Whitman_sims.results:
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    ax[0,0].plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uFastC']/totalC[0]*100, label = sim)
    #ax.plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uSlowC']/totalC[0]*100, color = 'blue', label = 'uSlowC')
    ax[1,0].plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uNecroC']/totalC[0]*100, label = sim)
    ax[0,1].plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uPyC']/totalC[0]*100, label = sim)
    ax[1,1].plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['uSlowC']/totalC[0]*100, label = sim)

#ax[0,0].set_xlabel('Time (days)')
ax[0,0].set_ylabel('Percent of total C')
ax[0,0].legend(loc='center left', prop={'size':6})
ax[0,0].set_title('uFastC', y=1, pad=10)
ax[1,0].set_xlabel('Time (days)')
ax[1,0].set_ylabel('Percent of total C')
#ax[1,0].legend(fontsize='small')
ax[1,0].set_title('uNecroC', y=1, pad=-10)
ax[0,1].set_xlabel('Time (days)')
#ax[0,1].set_ylabel('Percent of total C')
#ax[0,1].legend(fontsize='small')
ax[0,1].set_title('uPyC')
ax[1,1].set_xlabel('Time (days)')
#ax[1,1].set_ylabel('Percent of total C')
#ax[1,1].legend(fontsize='small')
ax[1,1].set_title('uSlowC', y=1, pad=-10)

pyplot.subplots_adjust(hspace = 0.2)

pyplot.show()




# Put the model and experimental results on the same figure.
fig,ax=pyplot.subplots(nrows=2,ncols=1,clear=True)
for sim in Whitman_sims.results:
    totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    if (sim=='no burn sandy soil' or sim=='high sev burn sandy soil'): 
        ax[0].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2']/totalC[0]*100, label = 'Model results: '+sim, linewidth=3.0) # Cumulative % of initial C respired
    else:
        ax[1].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2']/totalC[0]*100, label = 'Model results: '+sim, linewidth=3.0) # Cumulative % of initial C respired


# Add laboratory respiration data to figure:
no_burn_Pinus = df[(df["burn.trtmt"]=='control') & (df['vegetation']=='Pinus_banksiana')]
x_noburn_Pinus = no_burn_Pinus['whole.days.since.wet.up']
y_noburn_Pinus = no_burn_Pinus['new.cum_CO2C_g']/no_burn_Pinus['initial.total.C.g']*100
burn_30s_Pinus = df[(df["burn.trtmt"]=='low severity') & (df['vegetation']=='Pinus_banksiana')]
x_30burn_Pinus = burn_30s_Pinus['whole.days.since.wet.up']
y_30burn_Pinus = burn_30s_Pinus['new.cum_CO2C_g']/burn_30s_Pinus['initial.total.C.g']*100
burn_120s_Pinus = df[(df["burn.trtmt"]=='high severity') & (df['vegetation']=='Pinus_banksiana')]
x_120burn_Pinus = burn_120s_Pinus['whole.days.since.wet.up']
y_120burn_Pinus = burn_120s_Pinus['new.cum_CO2C_g']/burn_120s_Pinus['initial.total.C.g']*100
no_burn_Picea = df[(df["burn.trtmt"]=='control') & (df['vegetation']=='Picea_spp.')]
x_noburn_Picea = no_burn_Picea['whole.days.since.wet.up']
y_noburn_Picea = no_burn_Picea['new.cum_CO2C_g']/no_burn_Picea['initial.total.C.g']*100
burn_30s_Picea = df[(df["burn.trtmt"]=='low severity') & (df['vegetation']=='Picea_spp.')]
x_30burn_Picea = burn_30s_Picea['whole.days.since.wet.up']
y_30burn_Picea = burn_30s_Picea['new.cum_CO2C_g']/burn_30s_Picea['initial.total.C.g']*100
burn_120s_Picea = df[(df["burn.trtmt"]=='high severity') & (df['vegetation']=='Picea_spp.')]
x_120burn_Picea = burn_120s_Picea['whole.days.since.wet.up']
y_120burn_Picea = burn_120s_Picea['new.cum_CO2C_g']/burn_120s_Picea['initial.total.C.g']*100

ax[0].scatter(x_noburn_Pinus, y_noburn_Pinus, s=2, color = 'grey', label = "Lab results: No burn")
#ax[0].scatter(x_30burn_Pinus, y_30burn_Pinus, s=2, color = 'green', label = "Lab results: low severity burn")
ax[0].scatter(x_120burn_Pinus, y_120burn_Pinus, s=2, color = 'orange', label = "Lab results: high sev burn ")
ax[1].scatter(x_noburn_Picea, y_noburn_Picea, s=2, color = 'grey', label = "Lab results: No burn")
#ax[1].scatter(x_30burn_Picea, y_30burn_Picea, s=2, color = 'green', label = "Lab results: low severity burn")
ax[1].scatter(x_120burn_Picea, y_120burn_Picea, s=2, color = 'orange', label = "Lab results: high sev burn")

# Format axes and legend
ax[0].set_xlabel('Time (days)')
ax[0].set_ylabel('Cumulative C-CO$_2$ respired\n(% of initial C)')
ax[0].set_title('Sandy soil', y=1, pad=-15)
# ax[0].set_xlim([0,35])
#ax[0].legend(loc = 'upper left', prop={'size':8})
ax[1].set_xlabel('Time (days)')
ax[1].set_ylabel('Cumulative C-CO$_2$ respired\n(% of initial C)')
ax[1].set_title('Organic soil', y=1, pad=-15)
ax[1].legend(loc = 'upper left', prop={'size':8})
# ax[1].set_xlim([0,35])


pyplot.show() # Display plot
   
    
# Check model output by summing totalC and CO2
#fig,ax=pyplot.subplots(nrows=1)

# for sim in Whitman_sims.results:
#     totalC=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
#     ax.plot(Whitman_sims.t*365, Whitman_sims.results[sim][0]['CO2'] + totalC + Whitman_sims.results[sim][0]['livingMicrobeC'])

# ax.set_xlabel('Time (days)')
# ax.set_ylabel('Total C in model')
#ax.set_title('Laboratory incubation results', y=1, pad=-15)

#pyplot.show()


# Save results file
#df_output = Whitman_sims.results['burn'][0]
#df_output.to_csv('model_output/No_burn_Picea.csv', index_label='Time')

