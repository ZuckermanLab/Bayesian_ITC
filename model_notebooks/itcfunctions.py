import numpy as np
from scipy.optimize import root


#make list of injections
def get_injlist(injno,indiv_vol):
    injvol = [indiv_vol for i in range(injno)]
    return injvol

##function for root finding
def find_p_l(variables,pt,lt,k1,k2):
        p,l = variables
        
        equation_1 = pt - ((p) + ((2*l*p)/ k1 ) + (2*p*p*l) / (k1*k2))
        equation_2 = lt - ((l) + ((2*l*p)/ k1 ) + (p*p*l) / (k1*k2)) 

        return abs(equation_1),abs(equation_2)
        
##take volumes, initial concs (P0, Ls in bayesitc terms), return list of total concentrations
##based on the concentration calculation used in in Nguyen et. al. (2018).
def get_simulate_tots(V0,injvol,syringe_conc,cell_conc):

    injecting_tot = [0]
    cell_tot = [cell_conc]
    
    ##cumulative dilution factor
    dcum = 1  
    
    #calculate totals for each inj
    for inj in injvol:
        d  = 1 - (inj/V0)
        dcum *= d
        P = syringe_conc * (1-dcum)
        L = cell_conc * dcum
        
        ##change to numpy arrays when you get a chance
        injecting_tot.append(P)
        cell_tot.append(L)
    return injecting_tot,cell_tot



#function that takes in list of model parameters and calculates 
def get_dq_list(dg,ddg,dh,ddh,p_initial,l_initial,dh_0,inj_list,itc_constants):
    #unpack constants and convert dG to Kd
    Kb,T,V0 = itc_constants
    k1 = np.exp(dg/(Kb*T)) 
    k2 = np.exp((dg+ddg)/(Kb*T))
    
    #list of total concentrations 
    pt_list,lt_list = get_simulate_tots(V0,inj_list,p_initial,l_initial)
    q_list = []
    
    ##calculate Q for each injection
    for i in range(len(pt_list)):
        pt = pt_list[i]
        lt = lt_list[i]
        
        #root finding -- calculates [P] and [L] 
        #Uses default start point for first couple injections, then uses sol from previous round as start
        if i < 3:
            sol = root(find_p_l,method='lm',x0=(1e-8,1e-8),args=(pt,lt,k1,k2),options={'ftol':1e-20})
        else:
            sol = root(find_p_l,method='lm',x0=(sol.x),args=(pt,lt,k1,k2),options={'ftol':1e-20})
        pfree = sol.x[0]
        lfree = sol.x[1]
        
        #checks that the roots converge to ~0. value at if statement can be changed for stricter or more lenient checks
        #this often fails a couple times during the model's initial search, but corrects itself as walkers move to areas of parameter space that match the isotherm
        #when no root is found, Q calculations are incorrect, which leads to inflated likelihood
        root_check = pt - pfree - ((pfree*lfree*2)/k1) - ((2*pfree*pfree*lfree)/(k1*k2))
        if root_check > 1e-10:
            print('Could not find root')
            print(dg,ddg,dh,ddh,i)
            print(sol.x)
        
        
        #calculating model heat
        pl = 2*pfree*lfree / k1
        pl2 = lfree*pfree**2 / (k1*k2)
        q = V0 * (dh*pl + ((dh+(dh+ddh))*pl2))
        q_list.append(q)
        
    #generate list for dQ (i.e. change in heat per injection)
    dq_list = np.zeros(len(q_list)-1)
    for i in range(len(q_list)-1):
        dq = q_list[i+1] - q_list[i] + inj_list[i] / V0 * ((q_list[i+1] + q_list[i])/2)
        ##unit conversion from kcal to ucal (dh_0 in ucal already)
        dq_list[i] = dq*1e9 + dh_0
    return dq_list


#build synthetic isotherm using parameters listed in the function. 
def get_synthetic_itc(seed):
    
    np.random.seed(seed)

    ##itc constants
    ##these will change the isotherm being generated, so when working with synthetic data, they should match the constants in the notebook
    Kb = 0.001987
    T = 273.15 + 25
    V0 = 1.42e-3
    itc_constants = [Kb,T,V0]
    p_initial_stated = 500e-6
    l_initial_stated = 17e-6
    
    #kcal units -- these define thermo params for the synthetic model
    dg = -7
    dh = -10
    ddg = -1
    ddh = -1.5
    dh_0 = 0
    
    #ucal units
    sigma = 0.2
    
    theta_true = [dg,dh,ddg,ddh,p_initial_stated,l_initial_stated,sigma]

    #build injection list -- liter units, currently set to be 1 2 uL injection then 34 10 uL injections
    injection_count = 34
    injection_vol = 6e-6
    
    inj_list_l = [2e-6]
    for i in range(injection_count):
        inj_list_l.append(injection_vol)
    print(inj_list_l)
    
    #total concentrations per injection
    ptot,ltot = get_simulate_tots(V0,inj_list_l,p_initial_stated,l_initial_stated)
    #dQ values from get_dq_list function above
    true_dq = get_dq_list(dg,ddg,dh,ddh,p_initial_stated,l_initial_stated,dh_0,inj_list_l,itc_constants)
    ##add noise to isotherm around gaussian with width of sigma
    dq_obs = true_dq + np.random.normal(loc=0,scale=sigma,size=np.size(true_dq))
    
    return true_dq, dq_obs,theta_true,inj_list_l
