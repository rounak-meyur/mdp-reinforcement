# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 10:23:35 2019

@author: rounak
"""

#case = "kundur"
case = "ieee68"
node_file = case + "-bus.csv"
edge_file = case + "-line.csv"
gen_file = case + "-gen.csv"

# Import modules for setting up PSSE/python
import sys,os
import numpy as np
import cPickle as pkl
import matplotlib.pyplot as plt
workPath = os.getcwd()

# Get the different directory loaction into different variables
pathBase = workPath + "\\Base Case\\"
pathLib = workPath + "\\Libraries\\"
pathBin = workPath + "\\Binary Output\\"

# User defined libraries
sys.path.append(pathLib)
import pyDynamicslib
reload(pyDynamicslib)



buses = pyDynamicslib.Bus(pathBase+case+'-bus.csv')
lines = pyDynamicslib.Line(pathBase+case+'-line.csv')
tsfrs = pyDynamicslib.Tx(pathBase+case+'-tx.csv')
machs = pyDynamicslib.Mach(pathBase+case+'-gen.csv')



sol = pyDynamicslib.PowerFlow(buses,lines,tsfrs)
sol.ACPF(tol=0.5)
buses.update(sol)


model = pyDynamicslib.Machine(machs,buses)
Yred = pyDynamicslib.red_Ybus(buses,lines,tsfrs,machs)
ref = 0
Asys = model.Perturb(Yred)
w,V = pyDynamicslib.SlowCoherency(Asys,5)
gen,area = pyDynamicslib.Group(V)

print (np.column_stack((gen,area))

