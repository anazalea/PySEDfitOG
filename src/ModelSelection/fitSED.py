# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 19:55:43 2016

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
                
Top level module for model fitting
"""

from __future__ import print_function
from __future__ import division

import sys
import os

import math
import numpy as np
import param5
import fit_tools

#params = param5.SetFitSedParams('/Users/anneya/gzhk/forPySEDfit/test.param')
def main(params):
    '''
    Given a parameter FitSed parameter set, read data as specified, fit to 
    specified model set, with specified fitting method, record results as
    requested.
    '''
    
    # Quickly check that there's no output_overwrite conflict to save everyone's time
    if not params['output_overwrite'] and os.path.isfile(params['output_file']):
        raise IOError('output_overwrite is set to false but output_file exists.')
    
    # Read models
    modelFlux = np.genfromtxt(params['model_file'],usecols=(params['model_flux_columns']))
    if params['model_flux_unit']=='mag':
        modelFlux = 10.**(modelFlux/-2.5)
    modelParamInfo = np.array(params['model_param'])
    modelParams = np.genfromtxt(params['model_file'],usecols=(modelParamInfo[:,0].astype(int)))
    
    # Read data
    dataFlux = np.genfromtxt(params['data_file'],usecols=(params['data_flux_columns']))
    dDataFlux = np.genfromtxt(params['data_file'],usecols=(params['data_error_columns']))
    dataParamInfo = np.array(params['data_param'])
    dataParams = np.genfromtxt(params['data_file'],usecols=(dataParamInfo[:,0].astype(int)))
    
    if params['data_flux_unit']=='mag':
            dataFlux = 10.**(dataFlux/-2.5)
            #dDataFlux += magSoftening
            dDataFlux = dataFlux*dDataFlux/1.086
            
    # Identify fitting function
    fitFunc,fitParamInfo = fit_tools.FindFitFunction(params)
    
    # Make dictionary of additional method-dependent fit arguments
    fitArgs = fit_tools.GetFitArgs(params,modelParams)
    
    # If we want oldschool MC, generate
    if params['mcits']!=1 and params['oldschoolmc']==True:
        mcDataFlux = np.array([]).reshape((0,len(dataFlux[0])))
        mcDataParams = np.array([]).reshape((0,len(dataParams[0])))
        for i in range(len(dataFlux)):
            myMCDataFlux = np.array(params['mcits']*[list(dataFlux[i])])
            myMCDataParams = np.array(params['mcits']*[list(dataParams[i])])
            ns = np.random.randn(np.product(np.shape(myMCDataFlux))).reshape(np.shape(myMCDataFlux))
            fluxPerts = dDataFlux[i] * ns
            mcDataFlux = np.r_[mcDataFlux,fluxPerts+myMCDataFlux]
            mcDataParams = np.r_[mcDataParams,myMCDataParams]
        dataFlux = mcDataFlux
        dDataFlux = np.median(dDataFlux)*np.ones(np.shape(dataFlux))
        dataParams = mcDataParams
            
    # Fit 
    fit = fitFunc(dataFlux,dDataFlux,modelFlux,**fitArgs)
    outModelParams = modelParams[fit[:,0].astype(int)]
            
    # Scale model params where appropriate
    for i in range(len(modelParamInfo)):
        if modelParamInfo[i][3]=='True':
            outModelParams[:,i] = np.log10(fit[:,1]*(10**outModelParams[:,i]))
                    
    # Write output
    fit_tools.WriteOutput(dataParams,dataParamInfo,outModelParams,modelParamInfo[:,[0,1,2]],fit,fitParamInfo,params['output_file'])
    print(str(len(dataFlux))+' BBSEDs fit. Output saved to '+params['output_file'])
    
    
if __name__ == "__main__":
    pfile = sys.argv[1]
    args = sys.argv[2:]
    params = param5.SetFitSedParams(pfile,args)
    main(params)

