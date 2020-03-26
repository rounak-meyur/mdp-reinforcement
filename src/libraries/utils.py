# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 18:48:23 2020

Author: Rounak Meyur
Description: Utilities to set options and settings for simulations
"""
import sys
from recordtype import recordtype as rt
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix


inputcase = rt("case",field_names=['baseMVA','bus','branch','gen','order'])
extr = rt("ExtIndex",field_names=['bus','branch','gen'])
intr = rt("IntIndex",field_names=['bus','branch','gen'])
status = rt("Status",field_names=['on','off'])      # Named tuple for branch, gen, bus
tmp1 = rt("Temp",field_names=['e2i','i2e','status']) # Named tuple for gen and bus
tmp2 = rt("Temp",field_names=['status']) # Named tuple for branch
order = rt("Order",field_names=['extr','intr','bus','gen','branch','state'])

#%% Define constants
# define bus types
PQ      = 1
PV      = 2
REF     = 3
NONE    = 4

# define the indices
BUS_I       = 0    # bus number (1 to 29997)
BUS_TYPE    = 1    # bus type (1 - PQ bus, 2 - PV bus, 3 - reference bus, 4 - isolated bus)
PD          = 2    # Pd, real power demand (MW)
QD          = 3    # Qd, reactive power demand (MVAr)
GS          = 4    # Gs, shunt conductance (MW at V = 1.0 p.u.)
BS          = 5    # Bs, shunt susceptance (MVAr at V = 1.0 p.u.)
BUS_AREA    = 6    # area number, 1-100
VM          = 7    # Vm, voltage magnitude (p.u.)
VA          = 8    # Va, voltage angle (degrees)
BASE_KV     = 9    # baseKV, base voltage (kV)
ZONE        = 10   # zone, loss zone (1-999)
VMAX        = 11   # maxVm, maximum voltage magnitude (p.u.)      (not in PTI format)
VMIN        = 12   # minVm, minimum voltage magnitude (p.u.)      (not in PTI format)

# included in opf solution, not necessarily in input
# assume objective function has units, u
LAM_P       = 13   # Lagrange multiplier on real power mismatch (u/MW)
LAM_Q       = 14   # Lagrange multiplier on reactive power mismatch (u/MVAr)
MU_VMAX     = 15   # Kuhn-Tucker multiplier on upper voltage limit (u/p.u.)
MU_VMIN     = 16   # Kuhn-Tucker multiplier on lower voltage limit (u/p.u.)


# define the indices
GEN_BUS     = 0    # bus number
PG          = 1    # Pg, real power output (MW)
QG          = 2    # Qg, reactive power output (MVAr)
QMAX        = 3    # Qmax, maximum reactive power output at Pmin (MVAr)
QMIN        = 4    # Qmin, minimum reactive power output at Pmin (MVAr)
VG          = 5    # Vg, voltage magnitude setpoint (p.u.)
MBASE       = 6    # mBase, total MVA base of this machine, defaults to baseMVA
GEN_STATUS  = 7    # status, 1 - machine in service, 0 - machine out of service
PMAX        = 8    # Pmax, maximum real power output (MW)
PMIN        = 9    # Pmin, minimum real power output (MW)
PC1         = 10   # Pc1, lower real power output of PQ capability curve (MW)
PC2         = 11   # Pc2, upper real power output of PQ capability curve (MW)
QC1MIN      = 12   # Qc1min, minimum reactive power output at Pc1 (MVAr)
QC1MAX      = 13   # Qc1max, maximum reactive power output at Pc1 (MVAr)
QC2MIN      = 14   # Qc2min, minimum reactive power output at Pc2 (MVAr)
QC2MAX      = 15   # Qc2max, maximum reactive power output at Pc2 (MVAr)
RAMP_AGC    = 16   # ramp rate for load following/AGC (MW/min)
RAMP_10     = 17   # ramp rate for 10 minute reserves (MW)
RAMP_30     = 18   # ramp rate for 30 minute reserves (MW)
RAMP_Q      = 19   # ramp rate for reactive power (2 sec timescale) (MVAr/min)
APF         = 20   # area participation factor

# included in opf solution, not necessarily in input
# assume objective function has units, u
MU_PMAX     = 21   # Kuhn-Tucker multiplier on upper Pg limit (u/MW)
MU_PMIN     = 22   # Kuhn-Tucker multiplier on lower Pg limit (u/MW)
MU_QMAX     = 23   # Kuhn-Tucker multiplier on upper Qg limit (u/MVAr)
MU_QMIN     = 24   # Kuhn-Tucker multiplier on lower Qg limit (u/MVAr)


# define the indices
F_BUS       = 0    # f, from bus number
T_BUS       = 1    # t, to bus number
BR_R        = 2    # r, resistance (p.u.)
BR_X        = 3    # x, reactance (p.u.)
BR_B        = 4    # b, total line charging susceptance (p.u.)
RATE_A      = 5    # rateA, MVA rating A (long term rating)
RATE_B      = 6    # rateB, MVA rating B (short term rating)
RATE_C      = 7    # rateC, MVA rating C (emergency rating)
TAP         = 8    # ratio, transformer off nominal turns ratio
SHIFT       = 9    # angle, transformer phase shift angle (degrees)
BR_STATUS   = 10   # initial branch status, 1 - in service, 0 - out of service
ANGMIN      = 11   # minimum angle difference, angle(Vf) - angle(Vt) (degrees)
ANGMAX      = 12   # maximum angle difference, angle(Vf) - angle(Vt) (degrees)

# included in power flow solution, not necessarily in input
PF          = 13   # real power injected at "from" bus end (MW)       (not in PTI format)
QF          = 14   # reactive power injected at "from" bus end (MVAr) (not in PTI format)
PT          = 15   # real power injected at "to" bus end (MW)         (not in PTI format)
QT          = 16   # reactive power injected at "to" bus end (MVAr)   (not in PTI format)

# included in opf solution, not necessarily in input
# assume objective function has units, u
MU_SF       = 18   # Kuhn-Tucker multiplier on MVA limit at "from" bus (u/MVA)
MU_ST       = 19   # Kuhn-Tucker multiplier on MVA limit at "to" bus (u/MVA)
MU_ANGMIN   = 20   # Kuhn-Tucker multiplier lower angle difference limit (u/degree)
MU_ANGMAX   = 21   # Kuhn-Tucker multiplier upper angle difference limit (u/degree)


#%% Functions
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
    df_bus = pd.read_csv(pathdir+"bus"+casename+".csv").to_numpy()
    df_branch = pd.read_csv(pathdir+"branch"+casename+".csv").to_numpy()
    df_gen = pd.read_csv(pathdir+"gen"+casename+".csv").to_numpy()
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
            o = order(extr=extr(bus=mpc.bus,branch=mpc.branch,gen=mpc.gen),
                      intr=intr(bus=mpc.bus,branch=mpc.branch,gen=mpc.gen),
                      bus=tmp1(e2i={},i2e={},status=status(on=[],off=[])),
                      gen=tmp1(e2i={},i2e={},status=status(on=[],off=[])),
                      branch=tmp2(status(on=[],off=[])),state='e')
        else:
            o = mpc.order
        
        # save data matrices with external ordering
        o.extr.bus = mpc.bus
        o.extr.branch = mpc.branch
        o.extr.gen = mpc.gen
        
        # check that all buses have a valid BUS_TYPE
        bt = mpc.bus[:,BUS_TYPE]
        err = np.where((bt != PQ) & (bt != PV) & (bt != REF) & (bt != NONE))
        if np.shape(err)[1]!=0:
            print('ext2int: list of buses',err,'have invalid BUS_TYPE')
        
        # determine which buses, branches, gens are connected & in-service
        n2i = {int(k):i for i,k in enumerate(mpc.bus[:,BUS_I])}
        
        bs = (bt != NONE)                                           # bus status
        o.bus.status.on = np.where(bs)[0].tolist().copy()           # connected
        o.bus.status.off = np.where(bs == False)[0].tolist().copy() # isolated
        
        g2bi = [n2i[int(i)] for i in mpc.gen[:,GEN_BUS]]
        gs = np.logical_and(mpc.gen[:,GEN_STATUS],bs[g2bi])         # gen status
        o.gen.status.on = np.where(gs)[0].tolist().copy()           # on and connected
        o.gen.status.off = np.where(gs==False)[0].tolist().copy()   # off or isolated
        
        f2bi = [n2i[int(i)] for i in mpc.branch[:,F_BUS]]
        t2bi = [n2i[int(i)] for i in mpc.branch[:,T_BUS]]
        logic_bus = np.logical_and(bs[f2bi],bs[t2bi])
        brs = np.logical_and(mpc.branch[:,BR_STATUS],logic_bus)         # branch status
        o.branch.status.on = np.where(brs)[0].tolist().copy()           # on and connected
        o.branch.status.off = np.where(brs==False)[0].tolist().copy()   # off and disconnected
    
        # delete stuff that is "out"
        if len(o.bus.status.off)!=0:
            mpc.bus = np.delete(mpc.bus,o.bus.status.off,0)
        if len(o.branch.status.off)!=0:
            mpc.branch = np.delete(mpc.branch,o.branch.status.off,0)
        if len(o.gen.status.off)!=0:
            mpc.gen = np.delete(mpc.gen,o.gen.status.off,0)
        
        # apply consecutive bus numbering
        o.bus.i2e = {int(i):int(bus) for i,bus in enumerate(mpc.bus[:,BUS_I])}
        o.bus.e2i = {int(bus):int(i) for i,bus in enumerate(mpc.bus[:,BUS_I])}
        
        mpc.bus[:,BUS_I] = [o.bus.e2i[i] for i in mpc.bus[:,BUS_I]]
        mpc.gen[:,GEN_BUS] = [o.bus.e2i[i] for i in mpc.gen[:,GEN_BUS]]
        mpc.branch[:,F_BUS] = [o.bus.e2i[i] for i in mpc.branch[:,F_BUS]]
        mpc.branch[:,T_BUS] = [o.bus.e2i[i] for i in mpc.branch[:,T_BUS]]
    
        # reorder gens in order of increasing bus number
        genext = np.arange(0,mpc.gen.shape[0],step=1,dtype=int)
        genint = np.argsort(mpc.gen[:,GEN_BUS])
        o.gen.e2i = dict(zip(genext,genint))
        o.gen.i2e = dict(zip(genint,genext))
        mpc.gen = mpc.gen[genint,:]
        
        # save data matrices with internal ordering
        o.intr.bus = mpc.bus
        o.intr.branch = mpc.branch
        o.intr.gen = mpc.gen
    
        o.state = 'i'
        mpc.order = o
    return mpc

def int2ext(mpc):
    """
    Converts internal to external bus numbering.
    
    If the input is a single MATPOWER case struct, then it restores all buses, 
    generators and branches that were removed because of being isolated or off-line, 
    and reverts to the original generator ordering and original bus numbering. This 
    requires that the 'order' field created by EXT2INT be in place.
    
    Example:
        mpc = int2ext(mpc);
        
    See also EXT2INT
    Parameters
    ----------
    mpc : TYPE
        DESCRIPTION.

    Returns
    -------
    mpc : TYPE
        DESCRIPTION.

    """
    if mpc.order == None:
        print("Error!!! Order record type not present for conversion.")
        sys.exit(0)
    o = mpc.order
    if o.state == 'i':
        # save data matrices with internal ordering & restore originals
        o.intr.bus = mpc.bus
        o.intr.branch = mpc.branch
        o.intr.gen = mpc.gen
        mpc.bus = o.extr.bus
        mpc.branch = o.extr.branch
        mpc.gen = o.extr.gen

        # update data (in bus, branch, and gen only)
        mpc.bus[o.bus.status.on, :] = o.intr.bus
        mpc.branch[o.branch.status.on,:] = o.intr.branch
        mpc.gen[o.gen.status.on,:] = o.intr.gen[list(o.gen.i2e.values()),:]

        # # revert to original bus numbers
        int_bus = mpc.bus[o.bus.status.on, BUS_I]
        int_fbus = mpc.branch[o.branch.status.on, F_BUS]
        int_tbus = mpc.branch[o.branch.status.on, T_BUS]
        int_gbus = mpc.gen[o.gen.status.on, GEN_BUS]
        mpc.bus[o.bus.status.on, BUS_I] = np.array([o.bus.i2e[i] for i in int_bus])
        mpc.branch[o.branch.status.on, F_BUS] = np.array([o.bus.i2e[i] for i in int_fbus])
        mpc.branch[o.branch.status.on, T_BUS] = np.array([o.bus.i2e[i] for i in int_tbus])
        mpc.gen[o.gen.status.on, GEN_BUS] = np.array([o.bus.i2e[i] for i in int_gbus])
        
        o.state = 'e'
        mpc.order = o
    else:
        print("Case already in external ordering.")
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
    
    gen_bus = gen[:,GEN_BUS]
    gen_status = gen[:,GEN_STATUS]
    
    # element i, j is 1 if, generator index j at bus index i is ON
    Cg = csr_matrix((gen_status,(gen_bus, range(ng))), shape=(nb,ng))
    # number of generators at each bus that are ON
    bus_gen_status = Cg.dot(np.ones(shape=(ng, 1)))

    # form index lists for slack, PV, and PQ buses
    bt = bus[:,BUS_TYPE]
    ref = [i for i in range(nb) if bt[i] == REF and bus_gen_status[i]]
    pv  = [i for i in range(nb) if bt[i] == PV and bus_gen_status[i]]
    pq  = [i for i in range(nb) if bt[i] == PQ or not bus_gen_status[i]]

    # pick a new reference bus if for some reason there is none (may have been shut down)
    if len(ref)==0:
        ref = [pv.pop(0)]     # use the first PV bus and remove it from pv list
    
    return ref,pv,pq


def printpf(result,mpopt):
    """
    

    Parameters
    ----------
    result : TYPE
        DESCRIPTION.
    mpopt : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    text = ""
    text += "\n================================================================================"
    text += "\n|     Bus Data                                                                 |"
    text += "\n================================================================================"
    return text