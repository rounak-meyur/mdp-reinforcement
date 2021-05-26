# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:23:35 2019

Author: Rounak
"""

# Import modules
import sys,os
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
import pyModelinglib
from pyPowerFlowlib import PowerFlow

#%% Load case data
case = "ieee68"
buses = pyModelinglib.Bus(inppath+case+'-bus.csv')
lines = pyModelinglib.Line(inppath+case+'-line.csv')
tsfrs = pyModelinglib.Tx(inppath+case+'-tx.csv')
machs = pyModelinglib.Mach(inppath+case+'-gen.csv')


#%% Run ACPF
sol = PowerFlow(buses,lines,tsfrs)
sol.ACPF(flat_start=True,tol=0.05)
buses.update(sol)

sys.exit(0)

#%% Plot the results
def plot_pfresults(M,nx,ny,var):
    """
    """
    n,k = M.shape
    fig = plt.figure(figsize=(100,(ny*100.0)/nx))
    for i in range(n):
        ax = fig.add_subplot(ny,nx,i+1)
        ax.plot(M[i,:],color='crimson',marker='*',linestyle='dashed')
        ax.set_xlabel("Iterations")
        ax.set_ylabel(var)
    name = var+"-iterations"
    fig.savefig("{}{}.png".format(figpath,name),bbox_inches='tight')
    return

Vmag = sol.V_iter
theta = sol.theta_iter

plot_pfresults(Vmag,17,4,'voltage')
plot_pfresults(theta,17,4,'angle')


mismatch = sol.mismatch_iter
plot_pfresults(mismatch,17,8,'mismatch')





