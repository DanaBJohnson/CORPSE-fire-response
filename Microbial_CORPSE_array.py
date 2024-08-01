# This file holds the definition and key functions of the actual CORPSE model

# List of parameters the model needs
expected_params={'vmaxref_1': 'Relative maximum enzymatic decomp rates for first microbe pool (length 3)',
                 'vmaxref_2': 'Relative maximum enzymatic decomp rates for second microbe pool (length 3)',
                 'vmaxref_3': 'Relative maximum enzymatic decomp rates for third microbe pool (length 3)',
                 'vmaxref_4': 'Relative maximum enzymatic decomp rates for fourth microbe pool (length 3)',
        	'Ea':	'Activation energy (length 3)',
        	'kC':	'Michaelis-Menton parameter (length 3)',
        	'gas_diffusion_exp': 'Determines suppression of decomp at high soil moisture',
            'substrate_diffusion_exp':'Determines suppression of decomp at low soil moisture',
        	'minMicrobeC':	   'Minimum microbe biomass (fraction of total C)',
        	'Tmic': 'microbe lifetime at 20C (years)',
        	'et_1':  'Fraction of turnover from first microbial pool not converted to CO2',
        'et_2':  'Fraction of turnover from second microbial pool not converted to CO2',
        'et_3':  'Fraction of turnover from third microbial pool not converted to CO2',
        'et_4':  'Fraction of turnover from fourth microbial pool not converted to CO2',
        	'eup_1': 'Carbon uptake efficiency (length 3)',
        'eup_2':'Carbon uptake efficiency for second microbe pool (length3)',
        'eup_3':'Carbon uptake efficiency for third microbe pool (length3)',
        'eup_4':'Carbon uptake efficiency for fourth microbe pool (length3)',
        	'tProtected':	'Protected C turnover time (years)',
        	'protection_rate':'Protected carbon formation rate (year-1) (length 3)',
            'new_resp_units':True,
            }

# Names of the C types. Can edit this to change the number and name of pools in edited model simulations
chem_types = ['Fast','Slow','Necro', 'Py']

microbial_pools = ['1']

# Makes a list of the pools that should actually be in the model, including both protected and unprotected states
expected_pools = ['u'+t+'C' for t in chem_types]+\
                 ['p'+t+'C' for t in chem_types]+\
                 ['MBC_'+m for m in microbial_pools]+\
                 ['CO2','originalC']

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
             vmaxref_1=[2500,600,2000]; Relative maximum enzymatic decomp rates
             vmaxref_2=[2500,600,2000]; Relative maximum enzymatic decomp rates
             vmaxref_3=[2500,600,2000]; Relative maximum enzymatic decomp rates
             vmaxref_4=[2500,600,2000]; Relative maximum enzymatic decomp rates
             Ea=[37e3,54e3,50e3];     Activation energy
             kC=[0.01,0.01,0.01];     Michaelis-Menton parameter
             gas_diffusion_exp=2.5;   Determines suppression of decomp at high soil moisture
             minMicrobeC=1e-3;       Minimum microbe biomass (fraction of total C)
             Tmic=0.15;        microbe turnover rate
             et_1=0.5;           Fraction of turnover not converted to CO2
             et_2=0.5;
             et_3=0.5;
             et_4=0.5;
             eup_1=[0.6,0.05,0.6];  Carbon uptake efficiency - slow growing microbes
             eup_2=[0.4,0.03,0.4];  Carbon uptake efficiency - fast growing microbes
             eup_3=[0.4,0.03,0.4];  Carbon uptake efficiency - fast growing microbes
             eup_4=[0.4,0.03,0.4];  Carbon uptake efficiency - fast growing microbes
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

    # Calculate maximum potential C decomposition rate
    decomp=decompRate(SOM,T,theta,params)

    deadmic_C_production=0
    microbeTurnover=0
    # microbe turnover
    for m in microbial_pools:
        et=params['et_'+m]
        if m=='1':
            # Only calculate turnover if MBC pool is greater than 0
            microbeTurnover_1=where((SOM['MBC_1']>0), (SOM['MBC_1']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'],0);   # kg/m2/yr
            # Calculate cumulative microbial turnover
            microbeTurnover=microbeTurnover+microbeTurnover_1
            # Calculate cumulative necromass production
            deadmic_C_production=deadmic_C_production+microbeTurnover_1*et # actual fraction of microbe turnover
        elif m=='2':
            microbeTurnover_2=where((SOM['MBC_2']>0),(SOM['MBC_2']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'],0);   # kg/m2/yr
            microbeTurnover=microbeTurnover+microbeTurnover_2
            deadmic_C_production=deadmic_C_production+microbeTurnover_2*et
        elif m=='3':
            microbeTurnover_3=where((SOM['MBC_3']>0),(SOM['MBC_3']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'],0);   # kg/m2/yr
            microbeTurnover=microbeTurnover+microbeTurnover_3
            deadmic_C_production=deadmic_C_production+microbeTurnover_3*et
        elif m=='4':
            microbeTurnover_4=where((SOM['MBC_4']>0),(SOM['MBC_4']-params['minMicrobeC']*(sumCtypes(SOM,'u')))/params['Tmic'],0);   # kg/m2/yr
            microbeTurnover=microbeTurnover+microbeTurnover_4
            deadmic_C_production=deadmic_C_production+microbeTurnover_4*et


    if isinstance(microbeTurnover,float):
        microbeTurnover=max(0.0,microbeTurnover)
    else:
        microbeTurnover[microbeTurnover<0.0]=0.0

    # Calculate fraction of microbial turnover for each microbial pool conerted to CO2 (maintenance respiration)
    # Calculate cumulative maintenance respiration
    # Calculate cumulatve CO2 produced (maintenance respiration + decomposition of each C pool)
    # Calculate microbial growth for each microbial pool
    microbeGrowth_1=0
    microbeGrowth_2=0
    microbeGrowth_3=0
    microbeGrowth_4=0
    CO2prod=0
    maintenanceResp=0
    for m in microbial_pools:
        eup=params['eup_'+m]
        et=params['et_'+m]
        if m=='1':
            maintenanceResp_1=microbeTurnover_1*(1.0-et)
            maintenanceResp=maintenanceResp+maintenanceResp_1
            CO2prod=CO2prod+maintenanceResp_1     # CO2 production and cumulative CO2 produced by cohort
            for t in chem_types:
                CO2prod=CO2prod+decomp[m][t]*(1.0-eup[t])
                microbeGrowth_1=microbeGrowth_1+decomp[m][t]*eup[t]
        elif m=='2':
            maintenanceResp_2=microbeTurnover_2*(1.0-et)
            maintenanceResp=maintenanceResp+maintenanceResp_2
            CO2prod=CO2prod+maintenanceResp_2
            for t in chem_types:
                CO2prod=CO2prod+decomp[m][t]*(1.0-eup[t])
                microbeGrowth_2=microbeGrowth_2+decomp[m][t]*eup[t]
        elif m=='3':
            maintenanceResp_3=microbeTurnover_3*(1.0-et)
            maintenanceResp=maintenanceResp+maintenanceResp_3
            CO2prod=CO2prod+maintenanceResp_3
            for t in chem_types:
                CO2prod=CO2prod+decomp[m][t]*(1.0-eup[t])
                microbeGrowth_3=microbeGrowth_3+decomp[m][t]*eup[t]
        elif m=='4':
            maintenanceResp_4=microbeTurnover_4*(1.0-et)
            maintenanceResp=maintenanceResp+maintenanceResp_4
            CO2prod=CO2prod+maintenanceResp_4
            for t in chem_types:
                CO2prod=CO2prod+decomp[m][t]*(1.0-eup[t])
                microbeGrowth_4=microbeGrowth_4+decomp[m][t]*eup[t]
            


    # Update protected carbon
    protectedCturnover = dict([(t,SOM['p'+t+'C']/params['tProtected']) for t in chem_types])
    protectedCprod =     dict([(t,SOM['u'+t+'C']*params['protection_rate'][t]*claymod) for t in chem_types])
        # Some slow, protected C is being formed, but none of the pools are undergoing turnover > 0. Where is the protected C coming from??
        

    derivs=SOM.copy()
    for k in derivs.keys():
        derivs[k]=0.0
        
    for m in microbial_pools:
        if m=='1': derivs['MBC_1']=microbeGrowth_1-microbeTurnover_1; #derivs['CO2']=CO2prod_1
        elif m=='2': derivs['MBC_2']=microbeGrowth_2-microbeTurnover_2; #derivs['CO2']=derivs['CO2']+CO2prod_2
        elif m=='3': derivs['MBC_3']=microbeGrowth_3-microbeTurnover_3; #derivs['CO2']=derivs['CO2']+CO2prod_3
        elif m=='4': derivs['MBC_4']=microbeGrowth_4-microbeTurnover_4; #derivs['CO2']=derivs['CO2']+CO2prod_4
    derivs['CO2']=CO2prod

    total_decomp={}  
    for t in chem_types:
        total_decomp['u'+t+'C']=0
        total_decomp['p'+t+'C']=0
        for m in microbial_pools:
            total_decomp['u'+t+'C']=total_decomp['u'+t+'C']+decomp[m][t]
        derivs['u'+t+'C']=-total_decomp['u'+t+'C']+protectedCturnover[t]-protectedCprod[t]
        derivs['p'+t+'C']=-total_decomp['p'+t+'C']+protectedCprod[t]-protectedCturnover[t]

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

    for m in microbial_pools:
        if m=='1':
            vmax_1=Vmax('vmaxref_1',T,params)
        if m=='2':
            vmax_2=Vmax('vmaxref_2',T,params)
        if m=='3':
            vmax_3=Vmax('vmaxref_3',T,params)
        if m=='4':
            vmax_4=Vmax('vmaxref_4',T,params)
            

    # Decomposition rate of each C type
    for m in microbial_pools:
        if m=='1':
            decompRate_1={}
            decompRate={}
        if m=='2': 
            decompRate_2={}
        if m=='3':
            decompRate_3={}
        if m=='4':
            decompRate_4={}
            
    # Skip the decomposition calculation if there is no carbon or no microbe biomass (to avoid dividing by zero)
    dodecomp=(sumCtypes(SOM,'u')!=0.0)&(theta!=0.0)&(SOM['MBC_1']!=0.0)
    ##### Why does theta equal two values?
    
    for m in microbial_pools:    
        for t in chem_types:
            if dodecomp.any():
                if m=='1':
                    drate_1=where(dodecomp,vmax_1[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['MBC_1']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['MBC_1'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
                    decompRate_1[t]=drate_1
                    decompRate[m]=decompRate_1
                elif m=='2':
                    drate_2=where(dodecomp,vmax_2[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['MBC_2']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['MBC_2'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
                    decompRate_2[t]=drate_2
                    decompRate[m]=decompRate_2
                elif m=='3':
                    drate_3=where(dodecomp,vmax_3[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['MBC_3']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['MBC_3'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
                    decompRate_3[t]=drate_3
                    decompRate[m]=decompRate_3
                elif m=='4':
                    drate_4=where(dodecomp,vmax_4[t]*theta**params['substrate_diffusion_exp']*(SOM['u'+t+'C'])*SOM['MBC_4']/(sumCtypes(SOM,'u')*params['kC'][t]+SOM['MBC_4'])*(1.0-theta)**params['gas_diffusion_exp']/aerobic_max,0.0)
                    decompRate_4[t]=drate_4
                    decompRate[m]=decompRate_4

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
#  Doesn't include living MBC or CO2
def sumCtypes(SOM,prefix):
    out=SOM[prefix+chem_types[0]+'C']
    if len(chem_types)>1:
        for t in chem_types[1:]:
            out=out+SOM[prefix+t+'C']

    return out
