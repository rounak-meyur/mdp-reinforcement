# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:23:35 2019

@author: rounak
"""

#case = "kundur"
case = "ieee68"


# Import modules for setting up PSSE/python
import sys,os
import numpy as np
import matplotlib.pyplot as plt


# Get the different directory loaction into different variables
workpath = os.getcwd()
rootpath = os.path.dirname(workpath)

inppath = workpath + "/case-data/"
pathlib = rootpath + "/libs/"
figpath = workpath + "/figs/"

# User defined libraries
sys.path.append(pathlib)
import pyDynamicslib




buses = pyDynamicslib.Bus(inppath+case+'-bus.csv')
lines = pyDynamicslib.Line(inppath+case+'-line.csv')
tsfrs = pyDynamicslib.Tx(inppath+case+'-tx.csv')
machs = pyDynamicslib.Mach(inppath+case+'-gen.csv')



sol = pyDynamicslib.PowerFlow(buses,lines,tsfrs)
sol.ACPF(tol=0.5)
buses.update(sol)


model = pyDynamicslib.Machine(machs,buses)
Yred = pyDynamicslib.red_Ybus(buses,lines,tsfrs,machs)
ref = 0
Asys = model.Perturb(Yred)
w,V = pyDynamicslib.SlowCoherency(Asys,5)
gen,area = pyDynamicslib.Group(V)

print(np.column_stack((gen,area)))