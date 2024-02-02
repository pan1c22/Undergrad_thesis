#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 09:38:13 2022

@author: noahbungart
"""

import os
path1 = 'param_change/'
isExist = os.path.exists(path1)

if not isExist:
  
  # Create a new directory because it does not exist 
  os.makedirs(path1)
  print("The new directory is created!")