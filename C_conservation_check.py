import Microbial_sims
import Microbial_CORPSE_array as Microbial_CORPSE_array
from matplotlib import pyplot
import numpy as np

# Plot total C for each simulation
fig,ax=pyplot.subplots(nrows=4,ncols=1,clear=True,num='CORPSE results')
p=0

num_micro_pools=int(Microbial_sims.num_micro_pools)
for sim in Microbial_sims.results:
    totalC_substrate=Microbial_CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'u')+Microbial_CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'p')
    C_respired=Microbial_sims.results[sim][0]['CO2']
    totalMBC=0
    for m in np.arange(1,num_micro_pools+1,1):
        m=str(m) 
        totalMBC=Microbial_sims.results[sim][0]['MBC_'+m]+totalMBC
    
    totalC=totalC_substrate+C_respired+totalMBC
    
    ax[p].plot(Microbial_sims.t*365,totalC, color='purple')
    p=p+1

ax[0].set_title('Total C: no burn sandy soil', y=1, pad=-15)
ax[1].set_title('Total C: high sev burn sandy soil', y=1, pad=-15)
ax[2].set_title('Total C: no burn org soil', y=1, pad=-15)
ax[3].set_title('Total C: high sev burn org soil', y=1, pad=-15)

ax[p-1].set_xlabel('Time (days)')
ax[p-1].legend(fontsize='small')
pyplot.show()




        
        

# Plot C in each pool for the no burn sandy soil simulation
fig,ax=pyplot.subplots(nrows=4,ncols=1,clear=True,num='CORPSE results')
sim='no burn sandy soil'
ax[0].plot(Microbial_sims.t*365,Microbial_CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'u'), color='green')
ax[1].plot(Microbial_sims.t*365,Microbial_CORPSE_array.sumCtypes(Microbial_sims.results[sim][0], 'p'), color='green')

totalMBC=0

for m in np.arange(1,num_micro_pools+1,1):
    m=str(m) 
    totalMBC=Microbial_sims.results[sim][0]['MBC_'+m]+totalMBC
    
totalC=totalC_substrate+C_respired+totalMBC

ax[2].plot(Microbial_sims.t*365,totalMBC, color='green')
ax[3].plot(Microbial_sims.t*365,Microbial_sims.results[sim][0]['CO2'], color='green')

ax[0].set_title('Pool size for '+sim+' simulation', y=1, pad=-15)
ax[0].set_ylabel('Unprot. C')
ax[1].set_ylabel('Prot. C')
ax[2].set_ylabel('MBC')
ax[3].set_ylabel('CO2-C')

# ax[1].set_title('Protected C', y=1, pad=-15)
# ax[2].set_title('Total MBC', y=1, pad=-15)
# ax[3].set_title('CO2', y=1, pad=-15)
ax[3].set_xlabel('Time (days)')
pyplot.show()




# Plot C in each microbial pool for the no burn sandy soil simulation
if (Microbial_sims.num_micro_pools)>1:
    fig,ax=pyplot.subplots(nrows=num_micro_pools,ncols=1,clear=True,num='CORPSE results')
    p=0
    sim='no burn sandy soil'
    for m in np.arange(1,num_micro_pools+1,1):
        m=str(m) 
        ax[p].plot(Microbial_sims.t*365, Microbial_sims.results[sim][0]['MBC_'+m], color='purple')
        ax[p].set_ylabel('MBC '+m)
        p=p+1
    ax[0].set_title('Pool size for '+sim+' simulation', y=1, pad=-15)
    ax[p-1].set_xlabel('Time (days)')
    ax[p-1].legend(fontsize='small')
else:
    fig,ax=pyplot.subplots(nrows=1,ncols=1,clear=True,num='CORPSE results')
    p=0
    sim='high sev burn sandy soil'
    ax.plot(Microbial_sims.t*365, Microbial_sims.results[sim][0]['MBC_1'], color='purple')
    ax.set_ylabel('MBC 1')
    p=p+1 
    ax.set_title('Pool size for '+sim+' simulation', y=1, pad=-15)
    ax.set_xlabel('Time (days)')
    ax.legend(fontsize='small')

pyplot.show()







# C is conserved in single microbial pool model: 
import Whitman_sims
import CORPSE_array as CORPSE_array
from matplotlib import pyplot

# Plot total C for each simulation
fig,ax=pyplot.subplots(nrows=4,ncols=1,clear=True,num='CORPSE results')
p=0

for sim in Whitman_sims.results:
    totalC_substrate=CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u')+CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p')
    C_respired=Whitman_sims.results[sim][0]['CO2']
    totalMBC=Whitman_sims.results[sim][0]['livingMicrobeC']
    
    totalC=totalC_substrate+C_respired+totalMBC
    
    ax[p].plot(Whitman_sims.t*365,totalC, color='darkblue')
    p=p+1


ax[0].set_title('Total C: no burn sandy', y=1, pad=-15)
ax[1].set_title('Total C: high sev burn sandy', y=1, pad=-15)
ax[2].set_title('Total C: no burn org', y=1, pad=-15)
ax[3].set_title('Total C: high sev burn sandy', y=1, pad=-15)

ax[p-1].set_xlabel('Time (days)')
ax[p-1].legend(fontsize='small')
pyplot.show()



# # Plot C in each pool for the no burn sandy soil simulation
# fig,ax=pyplot.subplots(nrows=4,ncols=1,clear=True,num='CORPSE results')
# p=0

# Whitman_sims.results['no burn sandy soil']
# ax[0].plot(Whitman_sims.t*365,CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'u'), color='maroon')
# ax[1].plot(Whitman_sims.t*365,CORPSE_array.sumCtypes(Whitman_sims.results[sim][0], 'p'), color='maroon')

# totalMBC=0

# totalMBC=Whitman_sims.results[sim][0]['livingMicrobeC']+totalMBC
    
# totalC=totalC_substrate+C_respired+totalMBC

# ax[2].plot(Whitman_sims.t*365,totalMBC, color='maroon')
# ax[3].plot(Whitman_sims.t*365,Whitman_sims.results[sim][0]['CO2'], color='maroon')

# ax[0].set_title('Pool size for no burn sandy soil simulation', y=1, pad=-15)
# ax[0].set_ylabel('Unprot. C')
# ax[1].set_ylabel('Prot. C')
# ax[2].set_ylabel('MBC')
# ax[3].set_ylabel('CO2-C')
# ax[3].set_xlabel('Time (days)')
# pyplot.show()
