#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 15:14:24 2020

@author: capekp00
"""
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import differential_evolution
from scipy.optimize import dual_annealing

#read data
d = pd.read_csv('/mnt/580CBE2464C5F83D/pracovni/data_statistika/Data_z_clanku/Marstorp1999/Marstorp1999.csv', sep=',')
d
#define the simplest model
def Twopool (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2];  CO2=y[3]
    #define parameters
    v=pars[0];   k=pars[1]; f=pars[2]; m=pars[3]
    Y=pars[4]
    #Define fluxes
    ##glucose uptake
    uptake=v*Bs*G/(k+G)
    ##transfer function
    transfer=f*Br
    #death rate
    death=m*Bs
    #Define derivatives
    dGdt=-uptake
    dBrdt=uptake-transfer
    dBsdt=transfer*Y-death
    dCO2dt=transfer*(1-Y)
    
    return dGdt, dBrdt, dBsdt, dCO2dt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:5]
    #conversion factors
    conversions=pars[5:7]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    #r=pars_model[2]*(y[:, 1])*(1-pars_model[4])
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           y[:, 3].reshape(len(d.Time),1),#respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[6]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i, 0])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212cumul).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)/5).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[6]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i, 0])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
pars_opta=differential_evolution(obj_fun, [(1e-3, 10), (1e-2, 100), (1e-3, 10), (1e-8, 1), (0, 1), (0, 1),(0, 1)], 
                                          polish=True, maxiter=1000000)
print(pars_opta)
print(goodness(pars_opta.x))

#Dual annealing algorithm
pars_optb=dual_annealing(obj_fun, [(1e-3, 10), (1e-3, 100), (1e-3, 10), (1e-8, 1), (0, 1), (0, 1),(0, 1)])
print(pars_optb)
print(goodness(pars_optb.x))


np.savetxt('opt_parsa.csv', pars_opta.x.reshape(1,7), delimiter=',')
np.savetxt('opt_parsb.csv', pars_optb.x.reshape(1,7), delimiter=',')

#Two pool model with maintenance
def Twopool_main (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2];     CO2=y[3]
    #define parameters
    v=pars[0];   k=pars[1]; f=pars[2]; 
    m=pars[3];   #h=pars[4]; 
    Y=pars[4]
    #Define fluxes
    ##glucose uptake
    uptake=v*Bs*G/(k+G)
    ##transfer function
    transfer=f*Br
    #maintenance
    main=Bs*m
    #growth
    growth=max(0, (transfer-main)*Y)
    #respiration
    resp=max(0, (transfer-main)*(1-Y))+min(main, transfer)
    #death
    #death=m*Bs
        
    
    #Define derivatives
    dGdt=-uptake
    dBrdt=uptake-transfer
    dBsdt=growth-main
    dCO2dt=resp
    
    return dGdt, dBrdt, dBsdt, dCO2dt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:5]
    #conversion factors
    conversions=pars[5:7]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    #r=(pars_model[2]*y[:, 1])*(1-pars_model[5])+y[:, 2]*pars_model[4]
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           y[:, 3].reshape(len(d.Time),1),#respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[6]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i, 0])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool_main, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212cumul).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)/5).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[6]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i, 0])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool_main, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
pars_opt_maina=differential_evolution(obj_fun, [(1e-3, 10), (1e-2, 100), (1e-3, 10), (1e-8, 1), (0, 1), (0, 1), (0, 1)], 
                                          polish=True, maxiter=1000000)
print(pars_opt_maina)
print(goodness(pars_opt_maina.x))

#Dual annealing algorithm
pars_opt_mainb=dual_annealing(obj_fun, [(1e-3, 10), (1e-2, 100), (1e-3, 10), (1e-8, 1), (0, 1), (0, 1), (0, 1)])
print(pars_opt_mainb)
print(goodness(pars_opt_mainb.x))


np.savetxt('pars_opt_maina.csv', pars_opt_maina.x.reshape(1,7), delimiter=',')
np.savetxt('pars_opt_mainb.csv', pars_opt_mainb.x.reshape(1,7), delimiter=',')

#Two pool model with assimilation efficiency
def Twopool (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2]
    #define parameters
    v=pars[0];   k=pars[1]; f=pars[2]; 
    m=pars[3];   Ac=pars[4]; Y=pars[5]
    #Define fluxes
    ##glucose uptake
    uptake=v*Bs*G/(k+G)
    ##transfer function
    transfer=f*Br
    #death rate
    death=m*Bs
    #Define derivatives
    dGdt=-uptake
    dBrdt=Ac*uptake-transfer
    dBsdt=transfer*Y-death
    
    return dGdt, dBrdt, dBsdt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:6]
    #conversion factors
    conversions=pars[6:8]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    r=pars_model[2]*y[:, 1]*(1-pars_model[5])+y[:, 0]*y[:, 2]*pars_model[0]/(y[:, 0]+pars_model[1])*(1-pars_model[4])
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           r.reshape(len(d.Time),1),#respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[7]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA/5).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[7]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
pars_opt=differential_evolution(obj_fun, [(1e-3, 10), (1e-2, 100), (1e-3, 10), (1e-8, 1), (0, 1), (0, 1), (0, 1), (0, 1)], 
                                          polish=True, maxiter=1000000)
print(pars_opt)
print(goodness(pars_opt.x))

#Dual annealing algorithm
pars_opt=dual_annealing(obj_fun, [(1e-3, 10), (1e-2, 100), (1e-3, 10), (1e-8, 1), (0, 1), (0, 1), (0, 1), (0, 1)])
print(pars_opt)
print(goodness(pars_opt.x))


np.savetxt('opt_pars.csv', pars_opt.x.reshape(1,8), delimiter=',')

#Two pool model with linear uptake rate
def Twopool (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2]
    #define parameters
    v=pars[0];   f=pars[1]; 
    m=pars[2];   Y=pars[3]
    #Define fluxes
    ##glucose uptake
    uptake=v*Bs*G
    ##transfer function
    transfer=f*Br
    #death rate
    death=m*Bs
    #Define derivatives
    dGdt=-uptake
    dBrdt=uptake-transfer
    dBsdt=transfer*Y-death
    
    return dGdt, dBrdt, dBsdt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:4]
    #conversion factors
    conversions=pars[4:6]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    r=pars_model[1]*y[:, 1]*(1-pars_model[3])
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           r.reshape(len(d.Time),1),#respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[5]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)/5).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[5]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
pars_opt=differential_evolution(obj_fun, [(1e-3, 10), (1e-3, 10), (1e-8, 1), (0, 1), (0, 1), (0, 1)], 
                                          polish=True, maxiter=1000000)
print(pars_opt)
print(goodness(pars_opt.x))

#Dual annealing algorithm
pars_opt=dual_annealing(obj_fun, [(1e-3, 10),  (1e-3, 10), (1e-8, 1), (0, 1), (0, 1), (0, 1)])
print(pars_opt)
print(goodness(pars_opt.x))

np.savetxt('opt_pars.csv', pars_opt.x.reshape(1,6), delimiter=',')

#Two pool model with linear uptake rate and hyperbolic reserves release
def Twopool (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2]
    #define parameters
    v=pars[0];   f=pars[1];  kf=pars[2]
    m=pars[3];   Y=pars[4]
    #Define fluxes
    ##glucose uptake
    uptake=v*Bs*G
    ##transfer function
    transfer=f*Br/(kf+Br)
    #death rate
    death=m*Bs
    #Define derivatives
    dGdt=-uptake
    dBrdt=uptake-transfer
    dBsdt=transfer*Y-death
    
    return dGdt, dBrdt, dBsdt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:5]
    #conversion factors
    conversions=pars[5:7]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    r=pars_model[1]*y[:, 1]/(pars_model[2]+y[:, 1])*(1-pars_model[4])
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           r.reshape(len(d.Time),1),#respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[6]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA/5).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[6]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
pars_opt=differential_evolution(obj_fun, [(1e-3, 10), (1e-3, 10), (1e-3, 100), (1e-8, 1), (0, 1), (0, 1), (0, 1)], 
                                          polish=True, maxiter=1000000)
print(pars_opt)
print(goodness(pars_opt.x))

#Dual annealing algorithm
pars_opt=dual_annealing(obj_fun, [(1e-3, 10), (1e-3, 10), (1e-3, 100), (1e-8, 1), (0, 1), (0, 1), (0, 1)])
print(pars_opt)
print(goodness(pars_opt.x))

np.savetxt('opt_pars.csv', pars_opt.x.reshape(1,7), delimiter=',')

#Two pool model with hyperbolic uptake rate and reserves release
def Twopool (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2]
    #define parameters
    v=pars[0];   k=pars[1]; f=pars[2];  kf=pars[3]
    m=pars[4];   Y=pars[5]
    #Define fluxes
    ##glucose uptake
    uptake=v*Bs*G/(k+G)
    ##transfer function
    transfer=f*Br/(kf+Br)
    #death rate
    death=m*Bs
    #Define derivatives
    dGdt=-uptake
    dBrdt=uptake-transfer
    dBsdt=transfer*Y-death
    
    return dGdt, dBrdt, dBsdt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:6]
    #conversion factors
    conversions=pars[6:8]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    r=pars_model[2]*y[:, 1]/(pars_model[3]+y[:, 1])*(1-pars_model[5])
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           r.reshape(len(d.Time),1),#respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[7]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)/5).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[7]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
pars_opt=differential_evolution(obj_fun, [(1e-3, 10), (1e-3, 100), (1e-3, 10), (1e-3, 100), (1e-8, 1), (0, 1), (0, 1), (0, 1)], 
                                          polish=True, maxiter=1000000)
print(pars_opt)
print(goodness(pars_opt.x))

#Dual annealing algorithm
pars_opt=dual_annealing(obj_fun, [(1e-3, 10), (1e-3, 100), (1e-3, 10), (1e-3, 100), (1e-8, 1), (0, 1), (0, 1), (0, 1)])
print(pars_opt)
print(goodness(pars_opt.x))

np.savetxt('opt_pars.csv', pars_opt.x.reshape(1,8), delimiter=',')

#Two pool model with whatever
def Twopool (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2]
    #define parameters
    v=pars[0];   k=pars[1]; f=pars[2]; 
    m=pars[3];   a=pars[4]; Y=pars[5]
    #Define fluxes
    ##glucose uptake
    uptake=v*Bs*G/(k+G)
    ##transfer function
    transfer=f*Br
    #death rate
    death=m*Bs
    #Define derivatives
    dGdt=-uptake
    dBrdt=uptake-transfer+a*death
    dBsdt=transfer*Y-death
    
    return dGdt, dBrdt, dBsdt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:6]
    #conversion factors
    conversions=pars[6:8]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    r=pars_model[2]*y[:, 1]*(1-pars_model[5])
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           r.reshape(len(d.Time),1),#respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[7]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)/5).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[7]
    #Br_i = d.Cmicinit[0]/pars[5]     
    y0 = np.array([G_i, 0, Bs_i])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Twopool, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
pars_opt=differential_evolution(obj_fun, [(1e-3, 10), (1e-3, 100), (1e-3, 10), (1e-8, 10), (0, 1), (0, 0.8), (0, 1), (0, 1)], 
                                          polish=True, maxiter=1000000)
print(pars_opt)
print(goodness(pars_opt.x))

#Dual annealing algorithm
pars_opt=dual_annealing(obj_fun, [(1e-3, 10), (1e-3, 100), (1e-3, 10), (1e-8, 10), (0, 1), (0, 0.8), (0, 1), (0, 1)])
print(pars_opt)
print(goodness(pars_opt.x))

np.savetxt('opt_pars.csv', pars_opt.x.reshape(1,8), delimiter=',')


#Full soil biogeochemical model
def Full_model (y, t, pars):
    #define initial states
    G=y[0];    Br=y[1];    Bs=y[2];     DOC=y[3];       TOC=y[4];     ENZ=y[5];     CO2=y[6]
    #define fixed parameters
    vTOC=24 #maximal TOC decay rate [umol/umol(TOC)/d]
    kTOC=250/12.01*1000 #Affinity constant for TOC decay rate [umol/g]
    kDOC=0.26/12.01*1000 #Affinity constant for DOC uptake [umol/g]
    rENZ=0.024 #decay rate of enzymes [umol/umol(ENZ)/d]
    aBs=0.5    #proportion of Bs entering DOC
    Ye=0.75    #scaled yeild of extracellular enzymes production
    m=0.00028*24   #death rate 1/d
    aBr=0.75 #proportion of reserves used for the growth
    
    
    #define parameters
    v=pars[0]   #maximal substrate uptake rate [umol/umol(Bs)/d]
    kG=pars[1]  #Affinity constant for glucose uptake [umol/g]
    f=pars[2]   #mobilization rate of reserves  [umol/umol(Br)/d]
    #h=pars[3]   #maintenance coefficient [umol/umol(Bs)/d]
    Y=pars[3]   #Yield/Assimilation efficiency
    
    #Define fluxes
    ##TOC decay
    deTOC=vTOC*TOC*ENZ/(kTOC+TOC)
    ##ENZ decay
    deENZ=rENZ*ENZ
    ##DOC uptake
    uDOC=v*Bs*DOC/(kDOC+DOC)
    ##glucose uptake
    uG=v*Bs*G/(kG+G)
    ##transfer function
    transfer=f*Br
    ##growth
    growth=aBr*Y*transfer
    ##enzyme production
    eprod=(1-aBr)*Y*Ye*transfer
    ##death rate
    death=m*Bs
    ##respiration rate
    resp=aBr*(1-Y)*transfer+(1-aBr)*(1-Y*Ye)*transfer
    
    
    
    #Define derivatives
    dGdt=-uG
    dBrdt=uG+uDOC-transfer
    dBsdt=growth-death
    dDOCdt=deTOC+deENZ-uDOC+aBs*death
    dTOCdt=-deTOC+aBs*death
    dENZdt=eprod-deENZ
    dCO2dt=resp
    
    return dGdt, dBrdt, dBsdt, dDOCdt, dTOCdt, dENZdt, dCO2dt;
##use the model to calculate microbial biomass and respiration rate
def calc (model, pars, t, y0):
    #model parameters
    pars_model=pars[0:4]
    #conversion factors
    conversions=pars[4:6]
    #solve the model
    y=odeint(model,y0,t, args=(pars_model,))
    #calculate respiration rates
    #r=pars_model[2]*y[:, 1]*(1-pars_model[5])
    #calculate CFC
    fs=d.Cmicinit[0]*conversions[1]/d.DNAinit[0]
    CFC = conversions[0]*y[:, 1]+fs*y[:, 2]
    #calculate DNA
    DNA = conversions[1]*y[:, 2]
    #Create data with predictions
    yhat = np.concatenate((y[:, 0].reshape(len(d.Time),1),#glucose
                           y[:, 6].reshape(len(d.Time),1),#cumulative respiration
                           CFC.reshape(len(d.Time),1),
                           DNA.reshape(len(d.Time),1)), axis=1)
    return yhat

## Objective function is defined
def obj_fun (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[5]
    #Br_i = d.Cmicinit[0]/pars[5]
    DOC_i = d.Ctot[0]/12.01*1e4*0.001
    TOC_i = d.Ctot[0]/12.01*1e4
    ENZ_i = pars[6]
    y0 = np.array([G_i, 0, Bs_i, DOC_i, TOC_i, ENZ_i, 0])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Full_model, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    #weights
    weights=np.concatenate((np.nanmean(d.S).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.CO212cumul).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean((d.Cmic14+d.Cmic12)/5).repeat(len(d.Time)).reshape(len(d.Time),1),
                            np.nanmean(d.DNA).repeat(len(d.Time)).reshape(len(d.Time),1)), axis=1)
    out=np.nansum(((yhat_full-obs)/weights)**2)
    
    return out
## Goodness of fit
def goodness (x):
    #define parameters
    pars = x
    #initial conditions
    G_i = d.Sinit[0]
    Bs_i = d.DNAinit[0]/pars[5]
    #Br_i = d.Cmicinit[0]/pars[5]
    DOC_i = d.Ctot[0]/12.01*1e4*0.001
    TOC_i = d.Ctot[0]/12.01*1e4
    ENZ_i = pars[6]
    y0 = np.array([G_i, 0, Bs_i, DOC_i, TOC_i, ENZ_i, 0])
    #times
    t = d.Time
    #model simulations
    yhat_full = calc(Full_model, pars, t, y0)
    #observations
    obs=np.concatenate((np.array([d.S]).reshape(len(d.Time),1),
                        np.array([d.CO212cumul]).reshape(len(d.Time),1),
                        np.array([d.Cmic14+d.Cmic12]).reshape(len(d.Time),1),
                        np.array([d.DNA]).reshape(len(d.Time),1)), axis=1)
    R2=1-np.nansum((obs-yhat_full)**2)/np.nansum((obs-np.nanmean(obs))**2)
    ll=-np.nansum((obs-yhat_full)**2)/2/np.nanstd(obs)**2
    AIC = len(pars)*2 - 2*ll
    out = np.array([R2, ll, AIC])
    
    return out
## Parameters estimation
###Dual annealing algorithm
pars_optb=dual_annealing(obj_fun, [(1e-3, 10), #v
                                          (1e-2, 25), #kG
                                          (1e-3, 10), #f
                                          #(1e-8, 10), #h
                                          (0, 0.9),   #Y
                                          (0, 1), #fr
                                          (0, 1), #fd
                                          (0.0001, 0.42)], )#Enz0 
print(pars_optb)
print(goodness(pars_optb.x))

#np.savetxt('opt_pars_fulla.csv', pars_opt.x.reshape(1,9), delimiter=',')
np.savetxt('opt_pars_fullb.csv', pars_optb.x.reshape(1,7), delimiter=',')
