# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 17:44:08 2015

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
import init_PySEDfitOG
PySEDfitPath=init_PySEDfitOG.PySEDfitPath
import numpy as np
import spectrum
from astropy import units as u
from astropy import cosmology 
import dust_laws 
import filterv2
from copy import deepcopy 
import matplotlib.pyplot as plt
import IGM_attenuation


#------------------------------------------------------------------------------
def ProcessSpectrum(spec,dustLaw=None,ebvs=None,igmLaw=None,igmOpacities=None,cosmo=None,zs=None,filters=None,returnMag=True):
    '''
    Takes a single Spectrum instance and produces all combinations of reddening/redshifting requested, returns spectra or magnitudes
    Inputs:
        spec = base spectrum to be modified, instance of Spectrum
        dustLaw = name of dust law to be applied (or None), one of ['Calzetti2000','Calzetti1997','LMC','SMC','MW','Dor30']
        ebvs = list of E(B-V) values to be applied to base spectrum
        igmLaw = name of IGM attenuation law to be applied , one of ['madau',inoue']
        igmOpacities = list of IGM opacity values to be applied 
        cosmology = cosmology to use, instance of astropy.cosmology
        zs = list of redshifts to which base spectrum will be shifted
        filters = list of filters to be convolved with modified spectra, if None specified, full spectra are returned
        mags = bool, if True returns AB mags for each filter, if False returns erg s−1 cm−2 Hz−1
    Output:
        list of modified spectra, len(ebvs*igmOpacities*zs)
            if filters specified, list of bbsed objects at each combination of params
            if filters==None, full spectra are returned
    '''    
    # Validate input
    #print(spec.spec)
    
    ProcessSpectrumValidator(spec,zs,cosmology,ebvs,dustLaw,igmLaw,igmOpacities,filters)
    
    spectra = [] # spectra will store all of the modified spectra in a list
    if dustLaw == None:
        spectra = [deepcopy(spec)]

    # Add requested ISM reddening to base spectrum, each reddened spectrum is appended to spectra
    if dustLaw!=None:
        ebvs = np.asarray(ebvs)
        ebvs = ebvs.reshape((len(ebvs))) 
        for ebv in ebvs:
            newSpec = dust_laws.dustReddenSpectrum(spec,dustLaw,ebv)
            #print(newSpec.spec)
            spectra.append(newSpec)
    
    if igmLaw==None and cosmology==None: # Then we just want to return what we already have
        if filters==None: # Then we're returning spectra, don't need to do anything else
            return(spectra)
        else:
            bbseds = []
            for spect in spectra:
                mybbsed = spectrum.BBsed(params=deepcopy(spect.params))
                for filt in filters:
                    mag = spect.convolve(filt)
                    mybbsed.addMag(filt.name,mag)
                bbseds.append(deepcopy(mybbsed))
            return(bbseds)
            
    # If we got here, we either want to redshift and IGM attenuate or just redshift
    # We need a new list for storage (i.e. we don't necessarily want the z=0 cases from above)
    zSpectra = []
    bbseds = [] # which we'll fill if we're returning mags
    zs = np.asarray(zs)
    zs = zs.reshape(len(zs))
    
    if igmLaw!=None:
        igmOpacities = np.asarray(igmOpacities)
        igmOpacities = igmOpacities.reshape((len(igmOpacities)))
        
    for z in zs:
        dL = cosmo.luminosity_distance(z)
        for spect in spectra:
            mySpex = []
            zSpec = RedshiftSpectrum(spect,z=z,dL=dL,cosmo=cosmo)
            if igmLaw==None: # Or do we want to default to Madau?
                if filters==None:
                    zSpectra.append(zSpec)
                else:
                    mySpex.append(zSpec)
            elif igmLaw!=None:
                for igmO in igmOpacities:
                    igmSpec = IGM_attenuation.IGMAttenuateSpectrum(zSpec,igmLaw,igmO,z)
                    #igmSpec = deepcopy(zSpec) # will be igmSpec = IGM_attenuation.IGMattenuate()
                    # currently not available, but for each igmO we'll get a new spectrum which will either be appended to spectra or turned into a bbsed if we're not returning spectra
                    if filters==None:
                        zSpectra.append(igmSpec)
                    else:
                        mySpex.append(igmSpec)
            if filters!=None:
                for mySpec in mySpex:
                    mybbsed = spectrum.BBsed(params=deepcopy(mySpec.params))
                    for filt in filters:
                        filName =filt.name
                        mag = mySpec.convolve(filt,returnMag=returnMag)
                        mybbsed.addMag(filName,mag)
                    bbseds.append(deepcopy(mybbsed))
    if filters==None:
        return(zSpectra)
    else:
        return(bbseds)
#------------------------------------------------------------------------------    
def ProcessSpectrumValidator(spec,zs,cosmology,ebvs,dustLaw,igmLaw,igmOpacities,filters):
    # Spectrum
    try:
        if spec.type!='spectrum':
            raise TypeError('Input spectrum is not an instance of class Spectrum.')
    except: # if there's no type attribute of spec
        raise TypeError('Input spectrum is not an instance of class Spectrum.')
        
    # Redshifts
    zs = np.asarray(zs)
    zs = zs.reshape((len(zs)))
    if np.any(zs<0.):
        raise ValueError('All redshifts must be > 0.')
        
    try:
        h=cosmology.H0
    except:
        TypeError('Input cosmology must be an astropy.cosmology object.')
    
    # IGM
    if igmLaw!=None and zs==None:
        raise ValueError('If an IGM attenuation law is specified, redshifts must be specified.')
    if igmLaw!=None and igmLaw not in ['Madau','Inoue']:
        raise ValueError('The requested IGM attenuation law is not recognized.')
    # what constraints on IGM opacities? Just >0?
        
    # Dust
    if dustLaw!=None and dustLaw not in ['calzetti2000','calzetti1997','lmc','smc','mw','dor30']:
        raise ValueError('The requested dust law is not recognized.')
        ebvs = np.asarray(ebvs)
        ebvs = ebvs.reshape((len(ebvs)))
        if np.any(ebvs<0.):
            raise ValueError('Negative values of E(B-V) are not allowed.')
        if np.any(ebvs>1.0):
            raise ValueError('Values of E(B-V) > 1. are not allowed.')
            
    # Filters
    if filters!=None:
        if type(filters)!=list:
            raise TypeError('Filters must be specified in a list. If only one is specified, input as [filter].')
        for fil in filters:
            try:
                if fil.type!='filter':
                    print('All filters in filter list must be instances of class Filter.')
            except:
                print('All filters in filter list must be instances of class Filter.')
    
#------------------------------------------------------------------------------
def RedshiftSpectrum(spect,z,dL=None,cosmo=None):
    '''
    Redshift restframe spectrum to observed frame, apply cosmological dimming
    Inputs:
        spec = Rest frame spectrum to be redshifted, instance of class Spectrum
        z = redshift to which spectrum will be shifted
            one of [dL,cosmology] must be specified
        dL = luminosity distance associated with z, quantity (can be specified here to save time)
        cosmology = cosmology in which to calculate dL(z) if it hasn't been specified, astropy.cosmology instance
    Output:
        z_obs = Observed Spectrum, instance of class Spectrum
    '''
    if dL==None:
        dL = cosmo.luminosity_distance(z)
        print(dL)
    dL = dL.to('cm')

    #fNuNew = (spect.spec*u.Unit('m')**2)/(4*np.pi*dL**2/(1+z))
    uFnuN = u.Unit('erg')/u.Unit('s')/u.Unit('Hz')/u.Unit('cm')**2

    fNuNew = (spect.spec*(4*np.pi*((10 * u.Unit('parsec')).to('cm'))**2)/(4*np.pi*dL**2/(1+z))).to(uFnuN)

    newWavelengths = deepcopy(spect.wavelengths * (1.+z))
    newParams = deepcopy(spect.params)
    newParams['z']=z
    return(spectrum.Spectrum(newWavelengths.value,fNuNew.value,newWavelengths.unit,fNuNew.unit,params=newParams))
    
#------------------------------------------------------------------------------
#Test Stuff
'''
cosmo = cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

uf=np.genfromtxt('../../Testfiles/FTCs/uSDSS.ftc')
uSDSS = filterv2.Filter(uf[:,0],u.Unit('AA'),uf[:,1],'uSDSS')
gf=np.genfromtxt('../../Testfiles/FTCs/gSDSS.ftc')
gSDSS = filterv2.Filter(gf[:,0],u.Unit('AA'),gf[:,1],'gSDSS')
rf=np.genfromtxt('../../Testfiles/FTCs/rSDSS.ftc')
rSDSS = filterv2.Filter(rf[:,0],u.Unit('AA'),rf[:,1],'rSDSS')
iF=np.genfromtxt('../../Testfiles/FTCs/iSDSS.ftc')
iSDSS = filterv2.Filter(iF[:,0],u.Unit('AA'),iF[:,1],'iSDSS')
zf=np.genfromtxt('../../Testfiles/FTCs/zSDSS.ftc')
zSDSS = filterv2.Filter(zf[:,0],u.Unit('AA'),zf[:,1],'zSDSS')
filters = [uSDSS,gSDSS,rSDSS,iSDSS,zSDSS]
filNames = []
for filt in filters:
    filNames.append(filt.name)
lambdas = []
for filt in filters:
    lambdas.append(filt.effectiveWavelength().value)

f = np.genfromtxt('../../Testfiles/chab_tau0.1_m20.out.cols',usecols=(0,2))
spec = spectrum.Spectrum(f[:,0],f[:,1],u.Unit('Angstrom'),u.Unit('solLum')/u.Unit('Angstrom'))
spec.params['Age']='Yo mama.'
newspex=ProcessSpectrum(spec,cosmology=cosmo,zs=[1.0],igmLaw='Madau',igmOpacities=[0.25,0.5,1.0])

import matplotlib.cm as cm
i=0
for spex in newspex:
    plt.plot(spex.wavelengths.value,-2.5*np.log10(spex.spec.value)-48.6,color=cm.jet(i/len(newspex)),label='z='+str(round(spex.params['z'],3)))
    i+=1
plt.ylim(30,50)
plt.gca().invert_yaxis()
plt.xlim(500,9000)


newbbspex=ProcessSpectrum(spec,filters=filters,cosmology=cosmo,zs=[1.0],dustLaw='Calzetti2000',ebvs=[0.0,0.5])

i=0
for spex in newbbspex:
    pts = []
    for fn in filNames:
        pts.append(spex.mags[fn])
    plt.scatter(lambdas,pts,color=cm.jet(i/len(newspex)),s=10)
    i+=1
    
plt.legend()
'''