# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 21:43:08 2015

@author: anneya

                      .-. 
                 .--.(   ).--.                
      <-.  .-.-.(.->          )_  .--.
       `-`(     )-'             `)    )
         (o  o  )                `)`-'
        (      )                ,)   
        ( ()  )                 )
         `---"\    ,    ,    ,/`  
               `--' `--' `--'
                |  |   |   | 
                |  |   |   |
                '  |   '   |  
                
Functions to apply dust reddening curves to your sad blue :( spectra
"""

from __future__ import print_function
from __future__ import division

import sys
sys.path.append('.')
sys.path.append('../IO/')
import math
import numpy as np
import spectrum
from astropy import units as u
from copy import deepcopy 

def dustReddenSpectrum(spec,dustlaw,ebv):
    '''
    Apply dust reddening to rest frame galaxy spectrum
    
    Inputs:
        spec = spectrum object
        dustlaw = string, one of ['Calzetti2000','Calzetti1997','LMC','SMC','MW','Dor30']
        ebv = float in (0,1)
    Ouput:
        Astropy quantity of reddened spectrum of length len(spec.spec), with units of spec.spec
    '''
    #assert isinstance(spec,spectrum.Spectrum)
    assert isinstance(ebv,float)
    if ebv==0.0:
        newSpec = deepcopy(spec)
        newSpec.params['ebv']=0.0
        return(newSpec)
    dustLaws = ['calzetti2000','calzetti1997','lmc','smc','mw','dor30']
    dustLawFuncs = [Calzetti2000,Calzetti1997,LMC,SMC,MW,Dor30]
    if dustlaw not in dustLaws:
        raise ValueError(dustlaw,' is not a valid dust law.')
    if np.any(np.asarray(ebv)>=1.0):
        raise ValueError('Values of E(B-V) must be between 0 and 1.')
    nDustLaw = dustLaws.index(dustlaw)
    specUnit = spec.spec.unit
    wavelengthUnit = spec.wavelengths.unit
    reddenedSpec = dustLawFuncs[nDustLaw](1.0*spec.wavelengths.to('micron').value,spec.spec.value,ebv) * u.Unit(specUnit)  
    newSpec = spectrum.Spectrum(deepcopy(spec.wavelengths).value,reddenedSpec.value,wavelengthUnit,specUnit,params=deepcopy(spec.params)) 
    newSpec.params['ebv']=ebv
    return(newSpec)    
#------------------------------------------------------------------------------ 
def Calzetti2000(lam,lnu,ebv):
    '''
    Apply the Calzetti 20000 las as defined in Calzetti et al 2000 (ApJ 533 682) Equations 2,3,4
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''    
    r = 4.05
    k = 2.659 * (-1.857 + 1.040/lam) + r
    k[lam<0.63] = 2.659 * (-2.156 + (1.509/lam[lam<0.63]) - (0.198/lam[lam<0.63]**2) +\
                (0.011/lam[lam<0.63]**3)) + r
    newLnu = lnu * 10 ** (-0.4 * ebv * k)

    return(newLnu)
    
#------------------------------------------------------------------------------ 
def Calzetti1997(lam,lnu,ebv):
    '''
    Apply the Calzetti 1997 las as defined in Calzetti et al 1997 (astro-ph/9706121) Eqns 1 and 2. 
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''    
    k =  (((1.86 - 0.48/lam)/lam - 0.1) / lam) + 1.73
    k[lam<0.63] = 2.6596 * (-2.156 + (1.509/lnu[lam<0.63]) - (0.198/lnu[lam<0.63]**2) +\
                (0.011/lnu[lam<0.63]**3)) + 4.88
    newLnu = lnu * 10 ** (-0.4 * ebv * k)
    return(newLnu)
 
#------------------------------------------------------------------------------ 
def Fitzpatrick(lamInv,c1,c2,c3,c4,lambda0inv,gamma,ebv):
    '''
    Calculate tau values for one of the dust laws in the Fitzpatrick 1986 family
    Inputs:
            lamInv = NumPy array of wavelengths in microns^-1
            c1,c2,c3,c4 = Fitzpatricky constants (floats)
            lambda0inv = some other Fitzpatricky constant
            gamma = ditto that ^
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of exp(-1.*tau)
    '''   
    c4v = np.zeros(len(laminv))
    c4v[laminv>=5.9] = c4
    
    elvOverEbv = c1 + (c2*lamInv) + c3/((lamInv - lambda0inv**2/lamInv)**2 + gamma**2) +\
                c4v * (0.593*(lamInv-5.9)**2 + 0.0564*(lamInv-5.9)**3)
    avOverEbv = 3.1 # Assumes R(V)=A(V)/E(B-V)=3.1
    tau = (elvOverEbv + avOverEbv) * ebv/1.086
    return(np.exp(-1.*tau))
#------------------------------------------------------------------------------
def LMC(lam,lnu,ebv):
    '''
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)    
    Output:
            NumPy array of dust attenuated lnus
    '''   
    c1,c2,c3,c4 = -0.69,0.89,2.55,0.50
    lambda0inv = 4.608
    gamma = 0.994
    return(lam*Fitzpatrick(1./lam,c1,c2,c3,c4,lambda0inv,gamma,ebv))
#------------------------------------------------------------------------------  
def MW(lam,lnu,ebv):
    '''
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
    '''   
    c1,c2,c3,c4 = -0.38,0.74,3.96,0.26
    lambda0inv = 4.595
    gamma = 1.051
    return(lam*Fitzpatrick(1./lam,c1,c2,c3,c4,lambda0inv,gamma,ebv))
#------------------------------------------------------------------------------  
def Dor30(lam,lnu,ebv):
    '''
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
    Output:
            NumPy array of dust attenuated lnus
    '''   
    c1,c2,c3,c4 = -2.19,1.39,1.49,0.43
    lambda0inv = 4.606
    gamma = 0.894
    return(lam*Fitzpatrick(1./lam,c1,c2,c3,c4,lambda0inv,gamma,ebv))
#------------------------------------------------------------------------------ 
def SMC(lam,lnu,ebv):
    '''
    Apply the Prevot et al 1984 SMC law
    Inputs:
            lam = NumPy array of wavelengths in microns
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
    Output:
            NumPy array of dust attenuated lnus
    '''   
    lambdaInvPrevot = np.array([7.84, 7.52, 7.23, 6.98, 6.72,6.48, 6.27, 6.07, 5.88, \
    5.70, 5.52, 5.38, 5.24, 5.00, 4.73,4.50, 4.28, 4.09, 3.92, 3.75, 3.60, 3.46, \
    3.34, 3.22, 2.70,2.35, 1.89, 0])
    
    elvOverEbvPrevot = np.array([13.54, 12.52, 11.51, 10.80, 9.84, 9.28,9.06, 8.49, \
    8.01, 7.71, 7.17, 6.90, 6.76, 6.38, 5.85, 5.30,4.53, 4.24, 3.91, 3.49, 3.15, \
    3.00, 2.65, 2.29, 1.67, 1.00,0.00, 0.00])
    
    elvOverEbv = np.interp(1./lam,lambdaInvPrevot,elvOverEbvPrevot)
    avOverEbv = 3.1 # Assumes R(V)=A(V)/E(B-V)=3.1
    tau = (elvOverEbv + avOverEbv) * ebv/1.086
    return(lnu*np.exp(-1.*tau))
#------------------------------------------------------------------------------ 
    
    
    
    


