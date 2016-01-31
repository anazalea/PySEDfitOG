# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 18:59:59 2015

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
                
Model selection via chi2 using ball tree search when no MC iterations are requested
"""
from __future__ import print_function
from __future__ import division

import sys
sys.path.append('.')
import math
import numpy as np
import fit_tools 
from scipy.spatial.distance import seuclidean
from sklearn.neighbors import NearestNeighbors as nn
from sklearn.neighbors import BallTree

#------------------------------------------------------------------------------- 
def ChiTreeColor(dataFlux,dDataFlux,modelFlux):
    '''
    Finds the model that minimizes chi^2 distance for each object in color space 
    using a ball_tree search with seuclidean metric (weighting each dimension by 
    the variance in that flux)
    
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            
    Output:
            NumPy array of size (#objects,3)
            Columns [index of model with min chi^2, scale factor, chi^2 value]
    '''
    modelColors = modelFlux[:,1:] / modelFlux[:,:-1]
    dataColors = dataFlux[:,1:]/dataFlux[:,:-1]
    dDataColors = np.sqrt( (1./dataFlux[:,:-1])**2 * (dDataFlux[:,1:])**2 \
                + (dataFlux[:,1:]/dataFlux[:,:-1]**2)**2 * (dDataFlux[:,:-1])**2)
    results = np.array([]).reshape(0,3)
    for i in range(len(dataFlux)):
        tree = nn(n_neighbors=1,algorithm='ball_tree',metric='seuclidean',metric_params={'V':dDataColors[i]**2})
        tree.fit(modelColors)
        query = tree.kneighbors(dataColors[i],1)
        n,chi2 = query[1][0][0],query[0][0][0]**2.
        s = fit_tools.Scale(modelFlux[n],dataFlux[i],dDataFlux[i])
        results = np.r_[results,[[n,s,chi2]]]
    return(results)
#------------------------------------------------------------------------------ 
def ChiTreeFlux(dataFlux,dDataFlux,modelFlux):
    '''
    Finds the model that minimizes chi^2 distance for each object in flux space
    using a ball_tree search with seuclidean metric (weighting each dimension by 
    the variance in that flux)
    
    Inputs:
            dataFlux = observed fluxes, array of size (#objects,#filters)
            dDataFlux = flux uncertainties, array of size (#objects,#filters)
            modelFlux = fluxes of models, array of size (#models,#filters)
            
    Output:
            NumPy array of size (#objects,3)
            Columns [index of model with min chi^2, scale factor, chi^2 value]
    '''
    results = np.array([]).reshape(0,3)
    for i in range(len(dataFlux)):
        scales = fit_tools.Scale(modelFlux,dataFlux[i],dDataFlux[i])
        scaledModelFlux = (modelFlux.transpose() * scales.transpose()).transpose()
        tree = nn(n_neighbors=1,algorithm='ball_tree',metric='seuclidean',metric_params={'V':dDataFlux[i]**2})
        tree.fit(scaledModelFlux)
        query = tree.kneighbors(dataFlux[i],1)
        n,chi2 = query[1][0][0],query[0][0][0]**2.
        s = scales[int(n)]
        results = np.r_[results,[[n,s,chi2]]]
    return(results)
#------------------------------------------------------------------------------ 

