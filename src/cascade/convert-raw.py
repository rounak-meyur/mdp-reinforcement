# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:27:36 2020

@author: Rounak
"""

# Import modules
import sys,os
import pandapower.networks as pn
import pandapower as pp
import numpy as np
import matplotlib.pyplot as plt
print("Modules imported")

#%% Setup program

# Get the different directory loaction into different variables
workpath = os.getcwd()
rootpath = os.path.dirname(workpath)
inppath = rootpath + "/case-data/"
pathlib = rootpath + "/libs/"
figpath = workpath + "/figs/"

# User defined libraries
sys.path.append(pathlib)

#%% Load test case and run power flow
net = pn.case39()
df_buses = net.bus
df_lines = net.line
df_tsfr = net.trafo
df_gen = net.gen
df_load = net.load

pp.runpp(net,numba=False)

df_resgen = net.res_gen
df_resbus = net.res_bus
df_resline = net.res_line
df_restsfr = net.res_trafo
df_resload = net.res_load

#%% Fault
net.load.const_i_percent[0] = 100.0
net.load.p_mw *= 1000.0
df_load_new = net.load

pp.runpp(net,numba=False)
df_resgen_new = net.res_gen
df_resbus_new = net.res_bus
df_resline_new = net.res_line
df_restsfr_new = net.res_trafo
df_resload_new = net.res_load


























