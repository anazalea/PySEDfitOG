# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 10:42:51 2015

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

from __future__ import division, print_function
from astropy import constants as const
from astropy import units as u
import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os

class Filter:
    '''
    A filter.
    '''
    def __init__(self,x,xUnit,y,filtName):
        '''
        x = values of frequency or wavelength at which FTC is defined, array-like
        xUnit = astropy.units.Unit of type lengthor frequency
        y = filter transmission at x, array of length (len(x))
        '''
        self.type='filter'
        self.name=filtName
        wavelengths = np.asarray(x)
        ftc = np.asarray(y)
        if x.ndim!=1:
            raise ValueError('Wavelength/Frequency must be given as a 1d array.')
        if wavelengths.shape!=ftc.shape:
            raise ValueError('Wavelength/Frequency array and transmission array have different shapes.')
        if u.get_physical_type(xUnit) not in ['frequency','length']:
            raise TypeError(xUnit,' is neither a unit of frequency nor a unit of length. Get it together.')
        
        # Convert wavelengths/frequencies to microns
        wavelengths = x * xUnit
        wavelengths = wavelengths.to('Angstrom',equivalencies=u.spectral())
        
        # Eliminate any negative transmission values
        ftc = y
        ftc[y<0.]=0.
        # Reverse array orders if FTC was defined in frequency space
        if wavelengths[-1]<wavelengths[0]:
            wavelengths = np.flipud(wavelengths)
            ftc = np.flipud(ftc)
        self.wavelengths = wavelengths
        self.ftc = ftc
    
    def effectiveWavelength(self):
        return(self.wavelengths.unit*np.average(self.wavelengths,weights=self.ftc))

    