# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 22:14:36 2015

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
                
"""

from __future__ import print_function
from __future__ import division

import sys
sys.path.append('.')
import math
import numpy as np
from astropy import units as u
from copy import deepcopy 
import matplotlib.pyplot as plt
plt.rc('text',usetex=False)
plt.rc('font',family='serif')

dlaLS = np.genfromtxt('DLAcoeff.txt')
lafLS = np.genfromtxt('LAFcoeff.txt')
lambdaL = 911.8 # [AA]
    
#------------------------------------------------------------------------------ 
'''
def TauLAF_LS(lObs,zS,lJ,AJ1,AJ2,AJ3):
    
    #Calculate tau from Lyman Alpha Forest for line lJ
    #Inoue et al. 2014 eqn. X
    
    taus = 0.*lObs
    
    nonZero = (lJ<lObs)==(lObs<lJ*(1.+zS)) 
    
    whereAJ1 = (lObs < 2.2*lJ)==nonZero
    whereAJ2 = ((2.2*lJ <=lObs)==(lObs<5.7*lJ))==nonZero
    whereAJ3 = (5.7*lJ<=lObs)==nonZero
    
    taus[whereAJ1] = AJ1 * (lObs[whereAJ1]/lJ)**1.2
    taus[whereAJ2] = AJ2 * (lObs[whereAJ2]/lJ)**3.7
    taus[whereAJ3] = AJ3 * (lObs[whereAJ3]/lJ)**5.5

    return(taus)
'''
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
def tauInoue(lambdas,zS,nLines=39):
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
    
dlaLS = np.genfromtxt('DLAcoeff.txt')
lafLS = np.genfromtxt('LAFcoeff.txt')
lambdaL = 911.8 # [AA]

for z in np.arange(0,10.0,1.0):
    tauIGM(np.arange(0,100000.,10),z)
plt.xlim(0,13000)
plt.ylim(-0.1,1.1)
plt.legend(loc='lower right')
    
