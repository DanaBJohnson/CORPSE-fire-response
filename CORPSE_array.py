# This file holds the definition and key functions of the actual CORPSE model

# List of parameters the model needs
expected_params={	'vmaxref': 'Relative maximum enzymatic decomp rates (length 3)',
        	'Ea':	'Activation energy (length 3)',
        	'kC':	'Michaelis-Menton parameter (length 3)',
        	'gas_diffusion_exp': 'Determines suppression of decomp at high soil moisture',
            'substrate_diffusion_exp':'Determines suppression of decomp at low soil moisture',
        	'minMicrobeC':	   'Minimum microbial biomass (fraction of total C)',
        	'Tmic': 'Microbial lifetime at 20C (years)',
        	'et':		  'Fraction of turnover not converted to CO2',
        	'eup': 'Carbon uptake efficiency (length 3)',
        	'tProtected':	'Protected C turnover time (years)',
        	'protection_rate':'Protected carbon formation rate (year-1) (length 3)',
            'new_resp_units':True,
            }

# Names of the C types. Can edit this to change the number and name of pools in edited model simulations
chem_types = ['Fast','Slow','Necro', 'Py']

# Makes a list of the pools that should actually be in the model, including both protected and unprotected states
expected_pools = ['u'+t+'C' for t in chem_types]+\
                 ['p'+t+'C' for t in chem_types]+\
                 ['livingMicrobeC','CO2','originalC']

# Function for updating protected C formation rate for different soil textures, based on Mayes et al (2012) paper
#All soils: slope=0.4833,intercept=2.3282
#Alfisols: slope=0.5945, intercept=2.2788
def prot_clay(claypercent,slope=0.4833,intercept=2.3282,BD=1.15,porosity=0.4):
    ''' Calculate protection rate as a function of clay content, based on sorption isotherms from Mayes et al (2012) Table 3
    Calculates Qmax in mgC/kg soil from Mayes et al 2012, converted to g/m3 using bulk density
    Typically used as relative value for calculating protection_rate parameter.
    claypercent: Soil % clay (out of 100)
    slope: Either for all soils or a soil order, from Mayes paper
    intercept: Either for all soils or a soil order, from Mayes paper
    BD: Soil bulk density in g/cm3
    '''
    from numpy import log10
    prot=1.0*(10**(slope*log10(claypercent)+intercept)*BD*1e-6)
    return prot

# Check if the parameters sent to the model included the correct set of parameters and raise an error if not
def check_params(params):
    '''params: dictionary containing parameter values. Should contain these fields (showing reasonable default values):
             vmaxref=[2500,600,2000]; Relative maximum enzymatic decomp rates
             Ea=[37e3,54e3,50e3];     Activation energy
             kC=[0.01,0.01,0.01];     Michaelis-Menton parameter
             gas_diffusion_exp=2.5;   Determines suppression of decomp at high soil moisture
             minMicrobeC=1e-3;       Minimum microbial biomass (fraction of total C)
             Tmic=0.15;        Microbial turnover rate
             et=0.5;           Fraction of turnover not converted to CO2
             eup=[0.6,0.05,0.6];  Carbon uptake efficiency
             tProtected=75.0;     Protected C turnover time (years)
             protection_rate=[1.0,0.0,1.0];  Protected carbon formation rate (year-1)'''

    #from numpy import iterable,array
    unused_params=expected_params.copy()
    for k in params.keys():
        if k not in expected_params:
            raise ValueError('Parameter set contains unexpected parameter %s'%k)
        unused_params.pop(k)
        #if iterable(params[k]):
         #   params[k]=array(params[k])
    if len(unused_params)>0:
        for k in unused_params.keys():
            print ('Missing parameter: %s [%s]'%(k,unused_params[k]))
        raise ValueError('Missing parameters: %s'%unused_params.keys())

# The main model function. Given the current state of the model along with temperature, moisture, and parameters, it calculates the rate of change of all pools
from numpy import zeros,size,where,atleast_1d
def CORPSE_deriv(SOM,T,theta,params,claymod=1.0):
    '''Calculate rates of change for all CORPSE pools
       T: Temperature (K)
       theta: Soil water content (fraction of saturation)

       Returns same data structure as SOM'''

    theta=atleast_1d(theta)
    T=atleast_1d(T)

    theta[theta<0]=0.0
    theta[theta>1]=1.0

    et=params['et']
    eup=params['eup']

    # Calculate maximum potential C decomposition rate
    decomp=decompRate(SOM,T,theta,params)

    # Microbial turnover
    microbeTurnover=(SOM['livingMicrobeC']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'];   # kg/m2/yr
    if isinstance(microbeTurnover,float):
        microbeTurnover=max(0.0,microbeTurnover)
    else:
        microbeTurnover[microbeTurnover<0.0]=0.0

    maintenance_resp=microbeTurnover*(1.0-et)

    deadmic_C_production=microbeTurnover*et   # actual fraction of microbial turnover

    # CO2 production and cumulative CO2 produced by cohort
    CO2prod=maintenance_resp
    for t in chem_types:
        CO2prod=CO2prod+decomp[t]*(1.0-eup[t])

    microbeGrowth=CO2prod*0.0
    for t in chem_types:
        microbeGrowth=microbeGrowth+decomp[t]*eup[t]

    # Update protected carbon
    protectedCturnover = dict([(t,SOM['p'+t+'C']/params['tProtected']) for t in chem_types])
    protectedCprod =     dict([(t,SOM['u'+t+'C']*params['protection_rate'][t]*claymod) for t in chem_types])

    derivs=SOM.copy()
    for k in derivs.keys():
        derivs[k]=0.0
    derivs['livingMicrobeC']=microbeGrowth-microbeTurnover
    derivs['CO2']=CO2prod

    for t in chem_types:
        derivs['u'+t+'C']=-decomp[t]+protectedCturnover[t]-protectedCprod[t]
        derivs['p'+t+'C']=protectedCprod[t]-protectedCturnover[t]

    derivs['uNecroC']=derivs['uNecroC']+deadmic_C_production
    return derivs


# Decomposition rate
def decompRate(SOM,T,theta,params):

    # This only really needs to be calculated once
    # This corrects the units of vmaxref for the moisture function, so the 
    # units can be in actual 1/time. 
    # This essentially sets the maximum value of the moisture function to 1
    if params['new_resp_units']:
        theta_resp_max=params['substrate_diffusion_exp']/(params['gas_diffusion_exp']*(1.0+params['substrate_diffusion_exp']/params['gas_diffusion_exp']))
        aerobic_max=theta_resp_max**params['substrate_diffusion_exp']*(1.0-theta_resp_max)**params['gas_diffusion_exp']

    else:
        aerobic_max=1.0

    vmax=Vmax(T,params)

    # Decomposition rate of each C type
    decompRate={}
    # Skip the decomposition calculation if there is no carbon or no microbial biomass (to avoid dividing by zero)
    dodecomp=(sumCtypes(SOM,'u')!=0.0)&(theta!=0.0)&(SOM['livingMicrobeC']!=0.0)
    for t in chem_types:
        if dodecomp.any():
            drate=where(dodecomp,vmax[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['livingMicrobeC']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['livingMicrobeC'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
        decompRate[t]=drate

    return decompRate

def Vmax(T,params):
    '''Vmax function, normalized to Tref=293.15
    T is in K'''

    Tref=293.15;
    Rugas=8.314472;

    from numpy import exp

    Vmax=dict([(t,params['vmaxref'][t]*exp(-params['Ea'][t]*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))) for t in chem_types]);
    return Vmax

# Add together the C types. prefix is for specifying protected or unprotected (p or u)
#  Doesn't include living MBC
def sumCtypes(SOM,prefix):
    out=SOM[prefix+chem_types[0]+'C']
    if len(chem_types)>1:
        for t in chem_types[1:]:
            out=out+SOM[prefix+t+'C']

    return out
