# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 17:15:39 2020

@author: rouna
"""

import numpy as np
import pandas as pd

import pandapower.control as control
import pandapower.networks as nw
import pandapower.timeseries as timeseries
from pandapower.timeseries.data_sources.frame_data import DFData

# load a pandapower network
net = nw.mv_oberrhein(scenario='generation')
# number of time steps
n_ts = 95
# load your timeseries from a file (here csv file)
# df = pd.read_csv("sgen_timeseries.csv")
# or create a DataFrame with some random time series as an example
df = pd.DataFrame(np.random.normal(1., 0.1, size=(n_ts, len(net.sgen.index))),
                  index=list(range(n_ts)), columns=net.sgen.index) * net.sgen.p_mw.values
# create the data source from it
ds = DFData(df)

# initialising ConstControl controller to update values of the regenerative generators ("sgen" elements)
# the element_index specifies which elements to update (here all sgens in the net since net.sgen.index is passed)
# the controlled variable is "p_mw"
# the profile_name are the columns in the csv file (here this is also equal to the sgen indices 0-N )
const_sgen = control.ConstControl(net, element='sgen', element_index=net.sgen.index,
                                  variable='p_mw', data_source=ds, profile_name=net.sgen.index)

# do the same for loads
# df = pd.read_csv("load_timeseries.csv")
# create a DataFrame with some random time series as an example
df = pd.DataFrame(np.random.normal(1., 0.1, size=(n_ts, len(net.load.index))),
                  index=list(range(n_ts)), columns=net.load.index) * net.load.p_mw.values
ds = DFData(df)
const_load = control.ConstControl(net, element='load', element_index=net.load.index,
                                  variable='p_mw', data_source=ds, profile_name=net.load.index)

# starting the timeseries simulation for one day -> 96 15 min values.
timeseries.run_timeseries(net,numba=False)