# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 18:14:52 2020

Author: Rounak Meyur

Description: Models a Markov Decision Process for operating a power system. The csv
files are read first and the process is modeled with DC power flow assumptions.
"""

#%% Directories
import sys,os
workpath = os.getcwd()
basepath = workpath + "/case data/"
figpath = workpath + "/figures/"
libpath = workpath + "/libraries/"
sys.path.append(libpath)

from utils import loadcase
from powerutils import runpf
mpc = loadcase(basepath,'dat39')
mpc.bus[30,1] = 4

results = runpf(mpc)
