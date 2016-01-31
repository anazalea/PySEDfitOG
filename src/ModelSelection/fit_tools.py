# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 17:13:01 2015

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
                
Model selection tools common to all selection algorithms
"""

from __future__ import print_function
from __future__ import division

import sys
import math
import numpy as np

#-------------------------------------------------------------------------------   
def Scale(modelFlux,dataFlux,dDataFlux):
    '''
    Calculate scale factor of best fit according to Sawicki (2012) A5 
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of selected models, array of size (#objects,#filters)
            
    Output:
            Scale factors in numpy array of size (#objects)
            
    '''
    #  If there's only one object or one model, we need to reshape so 
    #  that array operations are working on the right indices
    if len(np.shape(modelFlux))==1:   
        modelFlux = np.reshape(modelFlux,(1,len(modelFlux)))
    if len(np.shape(dataFlux))==1:
        dataFlux = np.reshape(dataFlux,(1,len(dataFlux)))
    if len(np.shape(dDataFlux))==1:
        dDataFlux = np.reshape(dDataFlux,(1,len(dDataFlux)))
    s = np.sum(dataFlux*modelFlux/dDataFlux**2,axis=1) / np.sum(modelFlux**2/dDataFlux**2,axis=1)
    return(s)
   
#------------------------------------------------------------------------------ 
def FindFitFunction(params):
    '''
    Looks up the format of the output of the selected fitting method and gathers 
    parameters relevant to the fitting method from the control file in a dictionary to 
    be passed to the fitting function
    '''
    if params['fitting_method'] == 'brutedaisychain':
        import brute
        outParamInfo = np.array([['0','Best Fit Model Index','%i'],
                                 ['1','Scale Factor for Best fit Model','%f'],
                                 ['2','chi2','%f']])
        fitFunc = brute.BruteDaisyChain
       
    if params['fitting_method'] == "brutecolorspace":
        import brute
        outParamInfo = np.array([['0','Best Fit Model Index','%i'],
                                 ['1','Scale Factor for Best fit Model','%f'],
                                 ['2','chi2','%f']])
        fitFunc = brute.BruteColorSpace
        
    if params['fitting_method'] == 'brutefluxspace':
        import brute
        outParamInfo = np.array([['0','Best Fit Model Index','%i'],
                                 ['1','Scale Factor for Best fit Model','%f'],
                                 ['2','chi2','%f']])
        fitFunc = brute.BruteFluxSpace
    if params['fitting_method'] == 'brutefiterrorbars':
        import brute
        modelParamInfo = np.array(params['model_param'])
        outParamInfo = [['0','Best Fit Model Index','%i'],
                                 ['1','Scale Factor for Best fit Model','%f'],
                                 ['2','chi2','%f']] 
        i = 3
        for modelParam in modelParamInfo:
            newParam = [str(i),'Min. allowable '+modelParam[1]+' from dChi2',modelParam[2]]
            outParamInfo.append(newParam)
            i+=1
            newParam = [str(i),'Max. allowable '+modelParam[1]+' from dChi2',modelParam[2]]
            outParamInfo.append(newParam)
            i+=1
        fitFunc = brute.BruteFitErrorBars
    
    return([fitFunc,outParamInfo])
#------------------------------------------------------------------------------ 
def GetFitArgs(params,modelParams):
    if params['fitting_method'] == 'brutedaisychain':
        fitArgs = {}
       
    if params['fitting_method'] == "brutecolorspace":
        fitArgs = {}
        
    if params['fitting_method'] == 'brutefluxspace':
        fitArgs = {}
        
    if params['fitting_method'] == 'brutefiterrorbars':
        fitArgs = {'modelParams':modelParams,'dChi2':params['dchi2']}
    return(fitArgs)
       
    
    
#------------------------------------------------------------------------------  
def WriteHeader(dataParamInfo,modelParamInfo,fitParamInfo):
    '''
    Write an output file including each objects fit results and parameters of 
    the selected model
    Inputs:
            dataParams = 
            dataParamInfo = 
            outModelarams = 
            modelParamInfo = 
            fit = 
            fitParamInfo =
            Outfile = 
    Outputs:
            None
    '''
    head,fmt='',''
    i = 0
    for param in dataParamInfo:
        head+=' '+str(i)+' '+param[1]+'\n'
        fmt+=param[2]+' '
        i+=1
    for param in modelParamInfo:
        head+=' '+str(i)+' '+param[1]+'\n'
        fmt+=param[2]+' '
        i+=1
    for param in fitParamInfo:
        head+=' '+str(i)+' '+param[1]+'\n'
        fmt+=param[2]+' '
        i+=1

    return([head,fmt])

#------------------------------------------------------------------------------  
def WriteOutput(dataParams,dataParamInfo,outModelParams,modelParamInfo,fit,fitParamInfo,outfile):
    '''
    Write an output file including each objects fit results and parameters of 
    the selected model
    Inputs:
            dataParams = 
            dataParamInfo = 
            outModelarams = 
            modelParamInfo = 
            fit = 
            fitParamInfo =
            Outfile = 
    Outputs:
            None
    '''
    head,fmt = WriteHeader(dataParamInfo,modelParamInfo,fitParamInfo)
    outStuff = np.c_[dataParams,outModelParams,fit]

    np.savetxt(outfile,outStuff,fmt=fmt,header=head[:-1])
    return(None)
    
#------------------------------------------------------------------------------ 

