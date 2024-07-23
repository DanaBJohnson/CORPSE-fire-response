# This file holds the definition and key functions of the actual CORPSE model

# List of parameters the model needs
expected_params={'vmaxref_slow': 'Relative maximum enzymatic decomp rates for slower-growing microbial pool (length 3)',
                 'vmaxref_fast': 'Relative maximum enzymatic decomp rates for faster-growing microbial pool (length 3)',
        	'Ea':	'Activation energy (length 3)',
        	'kC':	'Michaelis-Menton parameter (length 3)',
        	'gas_diffusion_exp': 'Determines suppression of decomp at high soil moisture',
            'substrate_diffusion_exp':'Determines suppression of decomp at low soil moisture',
        	'minMicrobeC':	   'Minimum microbial biomass (fraction of total C)',
        	'Tmic': 'Microbial lifetime at 20C (years)',
        	'et_slow':  'Fraction of turnover not converted to CO2',
        'et_fast':  'Fraction of turnover from fast microbes not converted to CO2',
        	'eup_slow': 'Carbon uptake efficiency (length 3)',
        'eup_fast':'Carbon uptake efficiency for fast microbial pool (length3)',
        	'tProtected':	'Protected C turnover time (years)',
        	'protection_rate':'Protected carbon formation rate (year-1) (length 3)',
            'new_resp_units':True,
            }

# Names of the C types. Can edit this to change the number and name of pools in edited model simulations
chem_types = ['Fast','Slow','Necro', 'Py']

# Makes a list of the pools that should actually be in the model, including both protected and unprotected states
expected_pools = ['u'+t+'C' for t in chem_types]+\
                 ['p'+t+'C' for t in chem_types]+\
                 ['livingMicrobeC_fast','livingMicrobeC_slow','CO2','originalC']

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
             vmaxref_slow=[2500,600,2000]; Relative maximum enzymatic decomp rates
             vmaxref_fast=[2500,600,2000]; Relative maximum enzymatic decomp rates
             Ea=[37e3,54e3,50e3];     Activation energy
             kC=[0.01,0.01,0.01];     Michaelis-Menton parameter
             gas_diffusion_exp=2.5;   Determines suppression of decomp at high soil moisture
             minMicrobeC=1e-3;       Minimum microbial biomass (fraction of total C)
             Tmic=0.15;        Microbial turnover rate
             et=0.5;           Fraction of turnover not converted to CO2
             eup_slow=[0.6,0.05,0.6];  Carbon uptake efficiency - slow growing microbes
             eup_fast=[0.4,0.03,0.4];  Carbon uptake efficiency - fast growing microbes
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

    et_fast=params['et_fast']
    et_slow=params['et_slow']
    eup_slow=params['eup_slow']
    eup_fast=params['eup_fast']

    # Calculate maximum potential C decomposition rate
    decomp=decompRate(SOM,T,theta,params)

    # Microbial turnover
    microbeTurnover_slow=(SOM['livingMicrobeC_slow']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'];   # kg/m2/yr
    microbeTurnover_fast=(SOM['livingMicrobeC_fast']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'];   # kg/m2/yr
    microbeTurnover=(SOM['livingMicrobeC_fast']+SOM['livingMicrobeC_slow']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'];   # kg/m2/yr

    if isinstance(microbeTurnover,float):
        microbeTurnover=max(0.0,microbeTurnover)
    else:
        microbeTurnover[microbeTurnover<0.0]=0.0

    maintenance_resp_fast=microbeTurnover_fast*(1.0-et_fast)
    maintenance_resp_slow=microbeTurnover_slow*(1.0-et_slow)
    maintenance_resp=microbeTurnover_slow*(1.0-et_slow) + microbeTurnover_fast*(1.0-et_fast)

    deadmic_C_production=microbeTurnover_slow*et_slow + microbeTurnover_fast*et_fast   # actual fraction of microbial turnover

    # CO2 production and cumulative CO2 produced by cohort
    CO2prod_slow=maintenance_resp_slow
    CO2prod_fast=maintenance_resp_fast
    for t in chem_types:
        CO2prod_slow=CO2prod_slow+decomp['micro_slow'][t]*(1.0-eup_slow[t])
        CO2prod_fast=CO2prod_fast+decomp['micro_fast'][t]*(1.0-eup_fast[t])

    microbeGrowth_slow=CO2prod_slow*0.0
    microbeGrowth_fast=CO2prod_fast*0.0 # What is this doing? Creating empty variable    
    for t in chem_types:
        microbeGrowth_slow=microbeGrowth_slow+decomp['micro_slow'][t]*eup_slow[t]
        microbeGrowth_fast=microbeGrowth_fast+decomp['micro_fast'][t]*eup_fast[t]

    # Update protected carbon
    protectedCturnover = dict([(t,SOM['p'+t+'C']/params['tProtected']) for t in chem_types])
    protectedCprod =     dict([(t,SOM['u'+t+'C']*params['protection_rate'][t]*claymod) for t in chem_types])

    derivs=SOM.copy()
    for k in derivs.keys():
        derivs[k]=0.0
    derivs['livingMicrobeC_slow']=microbeGrowth_slow-microbeTurnover_slow
    derivs['livingMicrobeC_fast']=microbeGrowth_fast-microbeTurnover_fast
    derivs['CO2']=CO2prod_slow+CO2prod_fast

    for t in chem_types:
        derivs['u'+t+'C']=-decomp['micro_fast'][t]-decomp['micro_slow'][t]+protectedCturnover[t]-protectedCprod[t]
        
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

    vmax_slow=Vmax('vmaxref_slow', T,params)
    vmax_fast=Vmax('vmaxref_fast',T,params)

    # Decomposition rate of each C type
    decompRate_slow={}
    decompRate_fast={}
    decompRate={}
    # Skip the decomposition calculation if there is no carbon or no microbial biomass (to avoid dividing by zero)
    dodecomp=(sumCtypes(SOM,'u')!=0.0)&(theta!=0.0)&(SOM['livingMicrobeC_slow']!=0.0)
    for t in chem_types:
        if dodecomp.any():
            drate_slow=where(dodecomp,vmax_slow[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['livingMicrobeC_slow']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['livingMicrobeC_slow'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
            drate_fast=where(dodecomp,vmax_fast[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['livingMicrobeC_fast']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['livingMicrobeC_fast'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
        decompRate_slow[t]=drate_slow
        decompRate_fast[t]=drate_fast
        decompRate['micro_slow']=decompRate_slow
        decompRate['micro_fast']=decompRate_fast

    return decompRate

def Vmax(Micro_pool, T,params):
    '''Vmax function, normalized to Tref=293.15
    T is in K'''

    Tref=293.15;
    Rugas=8.314472;

    from numpy import exp

    Vmax=dict([(t,params[Micro_pool][t]*exp(-params['Ea'][t]*(1.0/(Rugas*T)-1.0/(Rugas*Tref)))) for t in chem_types]);

    return Vmax

# Add together the C types. prefix is for specifying protected or unprotected (p or u)
#  Doesn't include living MBC
def sumCtypes(SOM,prefix):
    out=SOM[prefix+chem_types[0]+'C']
    if len(chem_types)>1:
        for t in chem_types[1:]:
            out=out+SOM[prefix+t+'C']

    return out
