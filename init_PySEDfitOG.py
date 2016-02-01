# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 21:17:53 2016

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

import os
import sys
import init_PySEDfitOG
initpath = (init_PySEDfitOG.__file__).split('/')
global PySEDfitPath
PySEDfitPath = '/'.join(initpath[:-1])

sys.path.append(PySEDfitPath+'/src/SEDManipulation/')
sys.path.append(PySEDfitPath+'/src/IO/')

