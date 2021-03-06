# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 12:07:46 2015

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
print(PySEDfitPath)
import numpy as np
from copy import deepcopy 
import matplotlib.pyplot as plt
from astropy import units as u
from astropy import cosmology 
import filterv2
import spectrum
import param5
import redshift


def main(params):
    # Read spectra, depends on file type
    rfSpecs = ReadGalaxev(params['dotsed'],params['dotfourcolor'])
    
    # Read Filters
    filters = readFilters(params['filter_dir'],params['filter_names'])
    
    # Make cosmology
    print(params['cosmology'][0])
    if params['cosmology'][0]=='lcdm':
        cosmo = cosmology.LambdaCDM(Om0=params['cosmology'][1],Ode0=params['cosmology'][2],H0=params['cosmology'][3])
    else:
        wmaps,ns=[cosmology.WMAP5,cosmology.WMAP7,cosmology.WMAP9],[5.,7.,9.]
        cosmo = wmaps[ns.index(params['cosmology'][1])]

    # Parse zs
    if params['redshifts'][0]=='range':
        redshifts = np.arange(params['redshifts'][1],params['redshifts'][2],params['redshifts'][3])
    elif params['redshifts'][0]=='values':
        redshifts = np.array(params['redshifts'][1:])
    else:
        raise ValueError('Redshifts must be specified as a range or values.')
        
    # Dust Stuff
    if params['ebvs'][0]=='range':
        ebvs=np.arange(params['ebvs'][1],params['ebvs'][2],params['ebvs'][3])
    elif params['ebvs'][0]=='values':
        ebvs=params['ebvs'][1:]


    # make bbseds
    i=0
    for rfspec in rfSpecs:
        newbbspex = redshift.ProcessSpectrum(rfspec,filters=filters,cosmo=cosmo,\
            zs=redshifts,dustLaw=params['dust_law'],ebvs=ebvs,returnMag=params['returnmags'])
        if i==0:
            bbseds = newbbspex
            i+=1
        else:
            bbseds+=newbbspex

    # Write output header  
    head,sparams = WriteBBsedHeader(params,bbseds[0].params.keys(),params['returnmags'])
    bbsOut = []
    for bbs in bbseds:
        line=[]
        for sparam in sparams:
            line.append(bbs.params[sparam])
        for fil in params['filter_names']:
            line.append(bbs.mags[fil[0][:fil[0].index('.')]].value)
        bbsOut.append(line)
    bbsOut = np.array(bbsOut)
    np.savetxt(params['output_file'],bbsOut,header=head)
    return(head,bbsOut)
        
            
#------------------------------------------------------------------------------ 
def WriteBBsedHeader(params,specParams,returnMag):
    '''
    Writes header for .bbsed file where bbseds of processed spectra are recorder    
    Inputs:
        params = the input param list of makebbsed main(), ParamDict instance
        specParams = list of param names stored for each bbsed
        returnMag = bool, if True records units of 'mag' for each filter, if False, 'erg s−1 cm−2 Hz−1'
    Output:
        [header string, list of specParams(which will ne used to preserve order of params in output)]
    '''    
    head=' File generated by PySEDFit.makebbsed.py vX \n'
    head+=' full_param_list = '+str(params)+' \n'
    i=0
    if not returnMag:
        filtUnit = u.Unit('erg') * u.Unit('s')**-1 * u.Unit('cm')**-2 * u.Unit('Hz')**-1
    else:
        filtUnit = u.Unit('mag')
    for sp in specParams:
        head += ' '+str(i)+' '+sp+' \n'
        i+=1
    for filt in params['filter_names']:
        head += ' '+str(i)+ ' '+filt[0]+' '+str(filtUnit)+' \n'
        i+=1
    return([head[:-1],specParams])

#------------------------------------------------------------------------------ 
def readFilters(filterDir,filterStuff):
    '''
    Reads FTC information for each filter listed in param file, creates filter object for each
    returns list of filter objects
    Inputs:
        filterDir = location of FTCs. Each FTC file assumed to have 2 columns, 1:wavelength/frequency
                    2:transmission
        filterStuff = The list stored in the ['filter_names'] key of the paramDict created from
                    the .param file i.e. [[filtername,xUnitname],[filtername,xUnitname]...]
    Outputs:
        List of filter objects
    '''
    filters = []
    if filterDir[0]!='/':
        filterDir = '/'+filterDir
    if filterDir[-1]!='/':
        filterDir+='/'
    for fF in filterStuff:
        f = np.genfromtxt(filterDir+fF[0])
        try:
            filters.append(filterv2.Filter(f[:,0],u.Unit(fF[1]),f[:,1],fF[0][:fF[0].index('.')]))
        except:
            raise ValueError(fF)
    return(filters)

#------------------------------------------------------------------------------    
def ReadGalaxev(dotSed,dot4color):
    '''
    Reads data from a .sed and .4color file produced by galaxev (reformatted for SEDfit -WE SHOULD CHANGE THIS)
    Creates spectrum objects for each rest frame spectrum
    !!! DOES THIS RECORD THE PARAMS RIGHT IN ALL CASES?
    Inputs:
        dotSed = filename of .sed file
        dot4color = filename of .4color
        models = which models from galaxev files to be used on of:
                                ALL, [VALUES,1,2,3,...,n], [RANGE,min,max,step]
                                !!! SHOULD ADD SELECTION BY AGE AND STUFF HERE
    Output:
        List of spectrum objects
    '''    
    sed = np.genfromtxt(dotSed)
    fourColor = np.genfromtxt(dot4color)
    spex = []
    lSunA = u.Unit('Lsun')/u.Unit('AA')
    for i in range(len(fourColor)):
        ps = {'log(age/year)':fourColor[i][0],'SFR/yr':fourColor[i][-1],'Mstars':fourColor[i][6],'modelNumber':i}
        spex.append(spectrum.Spectrum(sed[:,0],sed[:,i+1],u.Unit('AA'),lSunA,params=ps))
    return(spex)
    


if __name__ == "__main__":
    pfile = sys.argv[1]
    args = sys.argv[2:]
    params = param5.SetMakeSedParams(pfile,args)
    main(params)

