# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:27:26 2015

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
import math
import numpy as np
from astropy import units as u
from copy import deepcopy 
import spectrum
import matplotlib.pyplot as plt

PySEDfitDir = '/Users/anneya/PySEDfit/'

def IGMAttenuateSpectrum(spec,igmLaw,igmOpacity,z):
    '''
    Apply IGM opacity to spectrum
    
    Inputs:
        spec = spectrum object
        dustlaw = string, one of ['Madau']
        igmOpacity = float > 0 
    Ouput:
        Astropy quantity of reddened spectrum of length len(spec.spec), with units of spec.spec
    '''
    #assert isinstance(spec,spectrum.Spectrum)
    assert isinstance(igmOpacity,float)
    igmLaws = ['Madau','Inoue']
    igmLawFuncs = [Madau,Inoue]
    if igmLaw not in igmLaws:
        raise ValueError(igmLaw,' is not a valid IGM Attenuation Law.')
    nLaw = igmLaws.index(igmLaw)
    specUnit = spec.spec.unit
    cosmicTrans = igmOpacity * igmLawFuncs[nLaw](spec.wavelengths.to('Angstrom').value,z)
    newSpec = spectrum.Spectrum(deepcopy(spec.wavelengths).value,cosmicTrans*deepcopy(spec.spec).value,spec.wavelengths.unit,u.Unit(specUnit),params=deepcopy(spec.params)) 
    newSpec.params['igmOpacity']=igmOpacity

    return(newSpec)    
#------------------------------------------------------------------------------ 
'''
def Madau(lam,lnu,z):
    
    Apply the Calzetti 20000 las as defined in Calzetti et al 2000 (ApJ 533 682) Equations 2,3,4
    Extrapolation below 1200A and above 2.2um
    Inputs:
            lam = NumPy array of wavelengths in A
            lnu = Numpy array of spectrum L_nu 
            ebv = E(B-V) a float in the range(0.,1.)
            
    Output:
            NumPy array of dust attenuated lnus
       
lam = spec.wavelengths.value
z=5.0
lymanLimit = 912 * u.Unit('AA')
lyLines = np.array([1216, 1026, 973, 950]).reshape((4,1))
Aj = np.array([0.0036, 0.0017, 0.0012, 0.00093]).reshape((4,1))

lamObs = lam * (1.+z)

# Line Blanketing
tauLines = Aj * (np.divide(np.vstack([lamObs,lamObs,lamObs,lamObs]),lyLines)) ** 3.46
mask = np.vstack([lam<lyLines[0][0],lam<lyLines[1][0],lam<lyLines[2][0],lam<lyLines[3][0]]).astype(float)
cosmicTrans = np.product(np.exp(-tauLines*mask),axis=0)
    
# Photoelectric
    
xEm = 1.+z
zC = lamObs / lymanLimit -1.
mask = lamObs > lymanLimit
'''
#########################################################
def Madau(lam,z):
    #lam = spec.wavelengths.value
    #z = 3.5
    wavEm = lam
    wavObs = lam * (1.+z)
    cosmicTrans = np.ones(len(lam))
    wavL = 912.
    aJ = [0.0036, 0.0017, 0.0012, 0.00093]
    wavJ = [1216., 1026., 973., 950.]
    
    for i in [0,1,2,3]:
        tauThisLine = aJ[i] * ((wavObs)/wavJ[i])**3.46
        mask = wavEm < wavJ[i]
        mask = wavObs < wavJ[i]*(1.+z)
        tauThisLine*=mask.astype(float)
        cosmicTrans*=np.exp(-1.*tauThisLine)
    
    
    xEm = 1. + z
    zC = wavObs/wavL - 1.
    xC = 1. + zC
    
    tauPhot = 0.25*xC**3*(xEm**0.46 - xC**0.46) \
            + 9.4*xC**1.5*(xEm**0.18 - xC**0.18) \
            - 0.7*xC**3*(xC**-1.32-xEm**-1.32) \
            - 0.023*(xEm**1.68 - xC**1.68)
            
    mask = wavObs<wavL*(1.+z)
    tauPhot *= mask.astype(float)
    #mask = wavObs<wavL
    #tauPhot *= mask.astype(float)
    cosmicTrans *= np.exp(-1.*tauPhot)
    #plt.plot(wavObs,cosmicTrans)
    return(cosmicTrans)


#########################################################
'''
Inoue 2014 Helpher functions
'''
dlaLS = np.genfromtxt(PySEDfitDir+'src/SEDManipulation/DLAcoeff.txt')
lafLS = np.genfromtxt(PySEDfitDir+'src/SEDManipulation/LAFcoeff.txt')
lambdaL = 911.8 # [AA]

#------------------------------------------------------------------------------ 

def TauLAF_LS(lObs,zS,lJ,AJ1,AJ2,AJ3):
  
    #Calculate tau from Lyman Alpha Forest for line lJ
    #Inoue et al. 2014 eqn. X
  
    taus = 0.*lObs
    
    nonZero = (lJ<lObs)==(lObs<lJ*(1.+zS))
    whereAJ1 = lObs < 2.2*lJ
    whereAJ2 = (2.2*lJ <=lObs)==(lObs<5.7*lJ)
    whereAJ3 = 5.7*lJ<=lObs
    
    taus[nonZero*whereAJ1.astype(bool)] = AJ1 * (lObs[nonZero*whereAJ1.astype(bool)]/lJ)**1.2
    taus[nonZero*whereAJ2.astype(bool)] = AJ2 * (lObs[nonZero*whereAJ2.astype(bool)]/lJ)**3.7
    taus[nonZero*whereAJ3.astype(bool)] = AJ3 * (lObs[nonZero*whereAJ3.astype(bool)]/lJ)**5.5

    return(taus)

#------------------------------------------------------------------------------         
def TauDLA_LS(lObs,zS,lJ,AJ1,AJ2):
    '''
    Calculate tau from Damped Lyman Alpha systems for line lJ
    Inoue et al. 2014 eqn. X
    '''
    taus = 0.*lObs
    
    nonZero = (lJ<lObs)==(lJ*(1.+zS))
    whereAJ1 = lObs<3.0*lJ
    whereAJ2 = lObs>=3.0*lJ
    
    taus[nonZero*whereAJ1.astype(bool)] = AJ1 * (lObs[nonZero*whereAJ1.astype(bool)]/lJ)**2.0
    taus[nonZero*whereAJ2.astype(bool)] = AJ2 * (lObs[nonZero*whereAJ2.astype(bool)]/lJ)**3.0
    
    return(taus)
#------------------------------------------------------------------------------
def TauLAF_LC(lObs,zS):
    '''
    
    '''
    lL = 911.8 #AA
    taus = 0.*lObs
    
    if zS < 1.2:
        p1 = lObs<lL*(1.+zS)
        taus[p1] = 0.325*( (lObs[p1]/lL)**1.2 - (1+zS)**-0.9 * (lObs[p1]/lL)**2.1)

    elif 1.2<=zS<4.7:
        p1 = lObs<2.2*lL
        taus[p1] = 2.55e-2 * (1+zS)**1.6 * (lObs[p1]/lL)**2.1 \
                        + 0.325*(lObs[p1]/lL)**1.2 - 0.250 * (lObs[p1]/lL)**2.1
                        
        p2 = (2.2*lL<=lObs)==(lObs<lL*(1.+zS))
        taus[p2] = 2.55e-2 * ((1+zS)**1.6 * (lObs[p2]/lL)**2.1 - (lObs[p2]/lL)**3.7 )
                        
    else: #  zS>=4.7: 
        p1 = lObs<2.2*lL
        taus[p1] = 5.22e-4 * (1+zS)**3.4*(lObs[p1]/lL)**2.1 + 0.325*(lObs[p1]/lL)**1.2 \
                                - 3.14e-2 * (lObs[p1]/lL)**2.1
                                
        p2 = (2.2*lL<=lObs)==(lObs<5.7*lL)
        taus[p2] = 5.22e-4 * (1+zS)**3.4 * (lObs[p2]/lL)**2.1 + 0.218*(lObs[p2]/lL)**2.1 \
                                            - 2.55e-2 * (lObs[p2]/lL)**3.7
                                            
        p3 = (5.7*lL<=lObs)==(lObs<lL*(1.+zS))
        taus[p3] = 5.22e-4 * ((1+zS)**3.4 *(lObs[p3]/lL)**2.1 - (lObs[p3]/lL)**5.5)
        
    return(taus)
#------------------------------------------------------------------------------
def TauDLA_LC(lObs,zS):
    '''
    
    '''
    lL=911.8 #AA
    taus = 0.*lObs
    
    if zS<2.0:
        p1 = lObs<lL*(1.+zS)
        taus[p1] = 0.211*(1.+zS)**2.0 -7.66e-2*(1.+zS)**2.3*(lObs[p1]/lL)**-0.3\
                                -0.135*(lObs[p1]/lL)**2.0
    else:
        p1 = lObs<3.0*lL
        taus[p1] = 0.634 + 4.70e-2 *(1+zS)**3.0 - 1.78e-2*(1+zS)**3.3 * (lObs[p1]/lL)**-0.3\
                            -0.135*(lObs[p1]/lL)**2.0 - 0.291*(lObs[p1]/lL)**-0.3
                            
        p2 = (3.0<=lObs)==(lObs<lL*(1.+zS))
        taus[p2] = 4.70e-2 * (1+zS)**3.0 - 1.78e-2*(1+zS)**3.3 \
                                * (lObs[p2]/lL)**-0.3 - 2.92e-2*(lObs[p2]/lL)**3.0     
    return(taus)
        
#------------------------------------------------------------------------------   
def Inoue(lambdas,zS,nLines=39):
    '''
    Calculates IGM attenuation at specified observed wavelengths for observations of
    a source at redshift zS according to Inoue et al. 2014
    Inputs:
            lambdas: NumPy array of wavelengths in Angstrom
            zS: float, redshift of course
            nLines: int, number of Lyman series lines to include (coefficients for 39 are included)
    '''
    # Lyman Series Contributions
    tauLAFls = 0.*lambdas
    tauDLAls = 0.*lambdas
    for j in range(nLines):
        tauLAFj = TauLAF_LS(lambdas,zS,lafLS[j][1],lafLS[j][2],lafLS[j][3],lafLS[j][4])
        tauLAFls += tauLAFj
        tauDLAj = TauDLA_LS(lambdas,zS,dlaLS[j][1],dlaLS[j][2],dlaLS[j][3])
        tauDLAls += tauDLAj
    
    # Lyman continuum
        
    tauLAFlc = TauLAF_LC(lambdas,zS)
    tauDLAlc = TauDLA_LC(lambdas,zS)
        
    tau = tauLAFls + tauDLAls + tauLAFlc + tauDLAlc
    plt.plot(lambdas,np.exp(-1.*tau),alpha=0.5,label='z='+str(zS))
    #plt.plot(lambdas,np.exp(-1.*tauLAFls),lw=7-z,alpha=0.25,label='z='+str(zS),color='black')
    #return([tauLAFls , tauDLAls , tauLAFlc , tauDLAlc])
    return(tau)   
