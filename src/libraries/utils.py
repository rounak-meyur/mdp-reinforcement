# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 18:48:23 2020

Author: Rounak Meyur
Description: Utilities to set options and settings for simulations
"""
import sys
from collections import namedtuple as nt
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix

inputcase = nt("case",field_names=['baseMVA','bus','branch','gen','order'])

ext = nt("ExtIndex",field_names=['bus','branch','gen'])
status = nt("Status",field_names=['on','off'])      # Named tuple for branch, gen, bus
tmp = nt("Temp",field_names=['e2i','i2e','status']) # Named tuple for gen and bus
order = nt("Order",field_names=['ext','bus','gen','branch','state'])



NONE = 4
REF = 3
PV = 2
PQ = 1

def loadcase(pathdir,casename = "dat39"):
    """
    Description
    ----------
    Returns a named tuple containing individual data frames as fields.

    Parameters
    ----------
    pathdir : TYPE string.
        DESCRIPTION path of directory containing the case file data.
    casename : TYPE string starting with 'dat', optional.
        DESCRIPTION. The default is "dat39" for the 39 bus system.

    Returns
    -------
    mpc : TYPE named tuple of type 'Case'.
        DESCRIPTION individual data frames are fields in the named tuple.
    """
    baseMVA = 100
    df_bus = pd.read_csv(pathdir+"bus"+casename+".csv")
    df_branch = pd.read_csv(pathdir+"branch"+casename+".csv")
    df_gen = pd.read_csv(pathdir+"gen"+casename+".csv")
    mpc = inputcase(baseMVA=baseMVA, bus=df_bus, branch=df_branch, gen=df_gen,
                    order=None)
    return mpc

def ext2int(mpc):
    """
    Generates the internal indexing of dataframes in the casedata named tuple.
    All isolated buses, off-line  generators and branches are removed along with any 
    generators or branches connected to isolated buses. Then the buses are renumbered 
    consecutively, beginning at 0, and the generators are sorted by increasing bus 
    number. Any 'ext2int' callback routines registered in the case are also invoked 
    automatically. All of the related indexing information and the original data 
    matrices are stored in an 'order' field in the struct to be used by INT2EXT to 
    perform the reverse conversions. If the case is already using internal numbering 
    it is returned unchanged.

    Parameters
    ----------
    mpc : TYPE named tuple of type WorkCase/InputCase
        DESCRIPTION case data: bus,branch,generator information along with order.

    Returns
    -------
    mpc: TYPE named tuple of type WorkCase/InputCase
        DESCRIPTION case data: bus,branch,generator information along with order.
    """
    
    if mpc.order == None or mpc.order.state=='e':
        if mpc.order == None:
            # initialize order named tuple
            o = order(ext=ext(bus=mpc.bus,branch=mpc.branch,gen=mpc.gen),
                      bus=tmp(e2i={},i2e={},status=status(on=[],off=[])),
                      gen=tmp(e2i={},i2e={},status=status(on=[],off=[])),
                      branch=status(on=[],off=[]),state='e')
        else:
            o = mpc.order
        
        # save data matrices with external ordering
        o.ext.bus.update(mpc.bus)
        o.ext.branch.update(mpc.branch)
        o.ext.gen.update(mpc.gen)
        
        # check that all buses have a valid BUS_TYPE
        bt = mpc.bus['BUS_TYPE'].values
        err = np.where((bt != PQ) & (bt != PV) & (bt != REF) & (bt == NONE))
        if np.shape(err)[1]!=0:
            print('ext2int: list of buses',err,'have invalid BUS_TYPE')
        
        # determine which buses, branches, gens are connected & in-service
        n2i = {k:i for i,k in enumerate(mpc.bus['BUS_I'].values)}
        
        bs = (bt != NONE)                               # bus status
        o.bus.status.on = np.where(bs)[0]               # connected
        o.bus.status.off = np.where(bs == False)[0]     # isolated
        
        g2bi = [n2i[i] for i in mpc.gen['GEN_BUS'].values]
        gs = mpc.gen['GEN_STATUS'] & bs[g2bi]           # gen status
        o.gen.status.on = np.where(gs)[0]               # on and connected
        o.gen.status.off = np.where(gs==False)[0]       # off or isolated
        
        f2bi = [n2i[i] for i in mpc.branch['F_BUS'].values]
        t2bi = [n2i[i] for i in mpc.branch['T_BUS'].values]
        brs = mpc.branch['BR_STATUS'] & bs[f2bi] & bs[t2bi] # branch status
        o.branch.status.on = np.where(brs)[0]               # on and connected
        o.branch.status.off = np.where(brs==False)[0]       # off and disconnected
    
        # delete stuff that is "out"
        if np.shape(o.bus.status.off)[1]!=0:
            mpc.bus.drop(mpc.bus.index[o.bus.status.off])
        if np.shape(o.branch.status.off)[1]!=0:
            mpc.branch.drop(mpc.branch.index[o.branch.status.off])
        if np.shape(o.gen.status.off)[1]!=0:
            mpc.gen.drop(mpc.gen.index[o.gen.status.off])
        
        # apply consecutive bus numbering
        o.bus.i2e = {i:bus for i,bus in enumerate(mpc.bus['BUS_I'].values)}
        o.bus.e2i = {bus:i for i,bus in enumerate(mpc.bus['BUS_I'].values)}
        
        mpc.bus['BUS_I'] = [o.bus.e2i[i] for i in mpc.bus['BUS_I']]
        mpc.gen['GEN_BUS'] = [o.bus.e2i[i] for i in mpc.gen['GEN_BUS']]
        mpc.branch['F_BUS'] = [o.bus.e2i[i] for i in mpc.branch['F_BUS']]
        mpc.branch['T_BUS'] = [o.bus.e2i[i] for i in mpc.branch['T_BUS']]
    
        # reorder gens in order of increasing bus number
        genext = mpc.gen['GEN_BUS'].values
        genint = np.argsort(genext)
        o.gen.e2i = dict(zip(genext,genint))
        o.gen.i2e = dict(zip(genint,genext))
        mpc.gen = mpc.gen.reindex(genint)
    
        o.state = 'i'
        mpc.order = o
    return mpc

def bustypes(bus,gen):
    """
    Builds index lists for each type of bus (REF, PV, PQ).
    [REF, PV, PQ] = BUSTYPES(BUS, GEN) 
    Generators with "out-of-service" status are treated as PQ buses with zero generation 
    (regardless of Pg/Qg values in gen). Expects BUS and GEN have been converted to 
    use internal consecutive bus numbering.

    Parameters
    ----------
    bus : TYPE pandas dataframe of bus data for the network
        DESCRIPTION Bus data with all fields.
    gen : TYPE pandas dataframe of generator data for the network
        DESCRIPTION Generator data with all fields.

    Returns
    -------
    ref : list of reference buses or slack buses.
    pv: list of PV buses
    pq: list of PQ buses

    """
    # get generator status
    nb = bus.shape[0]
    ng = gen.shape[0]
    
    gen_bus = gen['GEN_BUS'].values
    gen_status = gen['GEN_STATUS'].values
    
    # element i, j is 1 if, generator index j at bus index i is ON
    Cg = csr_matrix((gen_status,(gen_bus, range(ng))), shape=(nb,ng))
    # number of generators at each bus that are ON
    bus_gen_status = Cg.dot(np.ones(shape=(ng, 1)))

    # form index lists for slack, PV, and PQ buses
    ref = [i for i in range(nb) if bus['BUS_TYPE'][i] == REF and bus_gen_status[i]]
    pv  = [i for i in range(nb) if bus['BUS_TYPE'][i] == PV and bus_gen_status[i]]
    pq  = [i for i in range(nb) if bus['BUS_TYPE'][i] == PQ or not bus_gen_status[i]]

    # pick a new reference bus if for some reason there is none (may have been shut down)
    if len(ref)==0:
        ref = [pv.pop(0)]     # use the first PV bus and remove it from pv list
    
    return ref,pv,pq