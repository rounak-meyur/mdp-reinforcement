# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 18:48:23 2020

Author: Rounak Meyur
Description: Utilities to run power flow operations for power system analysis
"""
import sys,os
import time
import numpy as np
from recordtype import recordtype as rt
from scipy.sparse import csr_matrix, spdiags
from utils import loadcase,bustypes,ext2int,int2ext,printpf
from options import default_opt
import copy

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

#%% Recordtypes
zip_model = rt("ZIP",field_names=['z','i','p'])
solved_case = rt("solved", field_names=['mpc','success','iterations','time'])

#%% Functions

def dcpf(B,pbus,theta0,ref,pv,pq):
    """
    Description:
    -------
    dcpf solves the DC power flow problem for the case. The function solves for the 
    bus voltage angles at all but the reference bus.
    
    Parameters
    ----------
    B : numpy array or scipy sparse matrix with floating type entries
        the full system B matrix.
    pbus : numpy array of floating point entries
        vector of bus real power injections.
    theta0 : numpy array of floating point entries
        initial vector of bus voltage angles (in radians).
    ref : python list of integer entries
        lists of bus indices for the swing bus.
    pv : python list of integer entries
        lists of bus indices for the PV buses.
    pq : python list of integer entries
        lists of bus indices for the PQ buses.

    Returns
    -------
    theta: numpy array of floating point entries
        a vector of bus voltage angles (in radians).
    success: a boolean value
        1: Success
        0: Failure
    
    See also
    -------
    rundcpf, runpf
    """
    # arbitrary threshold on |Va| for declaring failure
    theta_threshold = 1e5

    # initialize result vector
    theta = theta0

    # update angles for non-reference buses
    Br_inv = np.linalg.inv(B[pv+pq,:][:,pv+pq].toarray())
    pinj = pbus[pv+pq] - (B[pv+pq,:][:,ref] * theta0[ref])
    theta[pv+pq] =  np.dot(Br_inv,pinj)

    # check for presence of *any* warning
    if max(abs(theta)) > theta_threshold:
        print("WARNING!!!: DC power flow failed to solve. Angles above threshold.")
        return theta0,False
    
    return theta,True


def rundcpf(casedata, mpopt, fname, solvedcase):
    """
    Description
    -------
    RESULTS = rundcpf(casedata, mpopt, fname, solvedcase)
    Runs a DC power flow, returning a RESULTS named tuple.
    
    Parameters
    ----------
    casedata : named tuple for the case data
        file with the case data.
    mpopt : named tuple for the options
        to override default options can be used to specify the solution algorithm, 
        output options termination tolerances.
    fname : string type data
        name of a file to which the pretty-printed output will be appended.
    solvedcase : name of file to which the solved case will be saved in MATPOWER case 
        format (M-file will be assumed unless the specified name ends with '.mat')

    Returns
    -------
    RESULTS : results named tuple, with the following fields:
        (all fields from the input MATPOWER case, i.e. bus, branch, gen, etc., 
         but with solved voltages, power flows, etc.) 
        order - info used in external <-> internal data conversion
        et - elapsed time in seconds
    
    See also
    -------
    runpf
    """
    # mpopt = mpoption(mpopt, 'model', 'DC') # Requires update for more options
    results = runpf(casedata, mpopt, fname, solvedcase)
    return results


def makeBdc(*args):
    """
    Builds the B matrices and phase shift injections for DC power flow.
    
    makeBdc(mpc)
    makeBdc(bus, branch)

    Returns the B matrices and phase shift injection vectors needed for a DC power 
    flow. The bus real power injections are related to bus voltage angles by
                P = BBUS * Va + PBUSINJ
    The real power flows at the from end the lines are related to the bus voltage 
    angles by
                Pf = BF * Va + PFINJ
    Does appropriate conversions to p.u. Bus numbers must be consecutive beginning at
    1 (i.e. internal ordering).
    
    Example:
        makeBdc(mpc)
        makeBdc(bus, branch)
       
    See also DCPF.


    Parameters
    ----------
    mpc : TYPE case data recordtype
        DESCRIPTION consists of bus and branch information.
        
    OR
    
    bus: TYPE Numpy 2-D array of bus data
        DESCRIPTION Data regarding buses in the power system. Should have internal
        ordering.
    branch: TYPE Numpy 2-D array of branch data
        DESCRIPTION Data regarding branches in the power sysem. Should have internal
        ordering.

    Returns
    -------
    Bbus: TYPE Scipy CSR sparse matrix of dimension nbus-by-nbus
        DESCRIPTION The bus susceptance matrix with tap changing ratios.
    Bf: TYPE Scipy CSR sparse matrix of dimension nbranch-by-nbus
        DESCRIPTION The weighted branch-bus incidence matrix with branch susceptance
        as the weights.
    Pbusinj: TYPE Numpy 1-D array of dimension nbus
        DESCRIPTION Bus injections due to phase shifting transformers.
    Pfinj: TYPE Numpy 1-D array of dimension nbranch
        DESCRIPTION power flow injections due to phase shifting transformers.

    """
    # extract from MPC if necessary
    if len(args) == 1:
        mpc     = args[0]
        bus     = mpc.bus
        branch  = mpc.branch
    elif len(args) == 2:
        bus = args[0]
        branch = args[1]
    else:
        print("Error!!! Wrong number of arguments passed.")
        sys.exit(0)

    # constants
    nb = bus.shape[0]          # number of buses
    nl = branch.shape[0]       # number of lines

    # check that bus numbers are equal to indices to bus (one set of bus numbers)
    if not np.array_equal(bus[:, BUS_I],np.array(range(nb))):
        print("makeBdc: buses must be numbered consecutively in bus matrix")
        print("\n Use ext2int() to convert to internal ordering")
        sys.exit(0)

    # for each branch, compute the elements of the branch B matrix and the phase
    # shift "quiescent" injections, where
    # 
    #       | Pf |   | Bff  Bft |   | Vaf |   | Pfinj |
    #       |    | = |          | * |     | + |       |
    #       | Pt |   | Btf  Btt |   | Vat |   | Ptinj |
    
    stat = branch[:, BR_STATUS]                     # ones at in-service branches
    b = np.divide(stat,branch[:, BR_X])             # series susceptance
    tap = np.ones(shape=(nl,))                      # default tap ratio = 1
    ind = np.where(branch[:, TAP])[0]               # indices of non-zero tap ratios
    tap[ind] = branch[ind, TAP]                     # assign non-zero tap ratios
    b = np.divide(b,tap)

    # build connection matrix Cft = Cf - Ct for line and from - to buses
    f = branch[:, F_BUS]           # list of "from" buses
    t = branch[:, T_BUS]           # list of "to" buses
    col_ind = np.concatenate((f,t))
    
    ind_br = np.array(range(nl))
    row_ind = np.concatenate((ind_br,ind_br))       # double set of row indices
    
    entry = np.ones(shape=(nl,))
    data = np.concatenate((entry,-entry))
    A = csr_matrix((data,(row_ind, col_ind)),shape=(nl,nb))       # incidence matrix
    
    # build Bf such that Bf * Va is the vector of real branch powers injected
    # at each branch's "from" bus
    bf_data = np.concatenate((b,-b))
    Bf = csr_matrix((bf_data,(row_ind, col_ind)),shape=(nl,nb))   # weighted incidence matrix
    
    # build Bbus
    Bbus = A.T * Bf
    
    # build phase shift injection vectors
    shift = branch[:, SHIFT] * np.pi/180
    Pfinj = np.multiply(b,-shift)
    Pbusinj = A.T * Pfinj
    
    return Bbus,Bf,Pbusinj,Pfinj


def makeSdzip(baseMVA, bus, mpopt=None):
    """
    Builds vectors of nominal complex bus power demands for ZIP loads.
    SD = MAKESDZIP(BASEMVA, BUS, MPOPT) returns a recordtype with three fields, each 
    an nb x 1 vectors. The fields 'z', 'i' and 'p' correspond to the nominal p.u. 
    complex power (at 1 p.u. voltage magnitude) of the constant impedance, constant 
    current, and constant power portions, respectively of the ZIP load model.
    
    Example:
        Sd = makeSdzip(baseMVA, bus, mpopt);


    Parameters
    ----------
    baseMVA : TYPE floating type data
        DESCRIPTION the base MVA for the system.
    bus : TYPE Numpy 2-D array
        DESCRIPTION Bus data information.
    mpopt : TYPE
        DESCRIPTION.

    Returns
    -------
    Sd: TYPE.

    """
    if mpopt!=None and mpopt.exp.sys_wide_zip_loads.pw != []:
        if len(mpopt.exp.sys_wide_zip_loads.pw) != 3:
            print("makeSdzip: 'exp.sys_wide_zip_loads.pw' must be a 3 element list")
        
        if abs(sum(mpopt.exp.sys_wide_zip_loads.pw) - 1) > 1e-10:
            print("makeSdzip: elements of 'exp.sys_wide_zip_loads.pw' must sum to 1");
        
        pw = mpopt.exp.sys_wide_zip_loads.pw
    else:
        pw = [1, 0, 0]
    
    if mpopt!=None and mpopt.exp.sys_wide_zip_loads.qw != []:
        if len(mpopt.exp.sys_wide_zip_loads.qw) != 3:
            print("makeSdzip: 'exp.sys_wide_zip_loads.qw' must be a 3 element list")
        
        if abs(sum(mpopt.exp.sys_wide_zip_loads.qw) - 1) > 1e-10:
            print("makeSdzip: elements of 'exp.sys_wide_zip_loads.qw' must sum to 1");
        
        qw = mpopt.exp.sys_wide_zip_loads.qw
    else:
        qw = pw
    
    Sd = zip_model(z=complex(0.0,0.0),i=complex(0.0,0.0),p=complex(0.0,0.0))
    Sd.z = (bus[:, PD] * pw[2]  + 1j * bus[:, QD] * qw[2]) / baseMVA
    Sd.i = (bus[:, PD] * pw[1]  + 1j * bus[:, QD] * qw[1]) / baseMVA
    Sd.p = (bus[:, PD] * pw[0]  + 1j * bus[:, QD] * qw[0]) / baseMVA
    return Sd


def makeSbus(*args):
    """
    Builds the vector of complex bus power injections.
    
    SBUS, DSBUS_DVM = MAKESBUS(BASEMVA, BUS, GEN)
    SBUS, DSBUS_DVM = MAKESBUS(BASEMVA, BUS, GEN, MPOPT, VM)
    SBUS, DSBUS_DVM = MAKESBUS(BASEMVA, BUS, GEN, MPOPT, VM, SG)
    
    returns the vector of complex bus power injections, that is, generation minus 
    load. Power is expressed in per unit. If the MPOPT and VM arguments are present 
    it evaluates any ZIP loads based on the provided voltage magnitude vector. If VM 
    is empty, it assumes nominal voltage. If SG is provided, it is a complex ng x 1 
    vector of generator power injections in p.u., and overrides the PG and QG columns 
    in GEN, using GEN only for connectivity information.
    
    
    It also computes the partial derivative of the bus injections with respect to 
    voltage magnitude. If VM is empty, it assumes no voltage dependence and returns 
    a sparse zero matrix.
    
    See also MAKEYBUS.

    Parameters
    ----------
    *args : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # default inputs
    if len(args) == 3:
        baseMVA = args[0]
        bus = args[1]
        gen = args[2]
        mpopt = None
        Vm = np.array([])
        Sg = np.array([])
    elif len(args) == 5:
        baseMVA = args[0]
        bus = args[1]
        gen = args[2]
        mpopt = args[3]
        Vm = args[4]
        Sg = np.array([])
    elif len(args) == 6:
        baseMVA = args[0]
        bus = args[1]
        gen = args[2]
        mpopt = args[3]
        Vm = args[4]
        Sg = args[5]
    else:
        print("Error!!! makeSbus: Invalid number of arguments passed. Check syntax.")
    
    nb = bus.shape[0]
    
    # get load parameters
    Sd = makeSdzip(baseMVA, bus, mpopt);
    
    if len(Vm)==0:
        dSbus_dVm = csr_matrix((nb, nb))
    else:
        dSbus_dVm = -spdiags(Sd.i + 2*np.multiply(Vm,Sd.z), 0, nb, nb)
    
    # compute per-bus generation in p.u.
    on_genind = np.where(gen[:,GEN_STATUS])[0]
    gbus = gen[on_genind,GEN_BUS]
    ngon = on_genind.shape[0]
    # connection matrix element i, j is 1 if gen on(j) at bus i is ON
    Cg = csr_matrix((np.ones(shape=(ngon,)),(gbus, range(ngon))), shape=(nb, ngon))  
    if len(Sg)!=0:
        Sbusg = Cg * Sg[on_genind]
    else:
        Sbusg = Cg * (gen[on_genind, PG] + 1j * gen[on_genind, QG]) / baseMVA

    # compute per-bus loads in p.u.
    if len(Vm)==0:
        Vm = np.ones(shape=(nb,))
    Sbusd = Sd.p + np.multiply(Sd.i,Vm) + np.multiply(Sd.z, np.square(Vm))

    # form net complex bus power injection vector
    # (power injected by generators + power injected by loads)
    Sbus = Sbusg - Sbusd
    return Sbus, dSbus_dVm



def runpf(*args):
    """
    

    Parameters
    ----------
    mpc   : TYPE
        DESCRIPTION.
    mpopt : TYPE
        DESCRIPTION.
    fname : TYPE
        DESCRIPTION.
    
    Returns
    -------
    results: TYPE
        DESCRIPTION.

    """
    # read data
    if len(args)==1:
        mpc = copy.deepcopy(args[0])
        mpopt = default_opt
        fname = None
    elif len(args)==2:
        mpc = copy.deepcopy(args[0])
        mpopt = args[1]
        fname = None
    elif len(args)==3:
        mpc = copy.deepcopy(args[0])
        mpopt = args[1]
        fname = args[2]
    else:
        print("Error!!! 'runpf': Invalid number of arguments.")
        sys.exit(0)
    
    # Create the results structure
    results = solved_case(mpc=mpc,success=False,iterations=0,time=0.0)
    
    # get bus index lists of each type of bus
    mpc = ext2int(mpc)
    baseMVA = mpc.baseMVA
    gen = mpc.gen
    bus = mpc.bus
    branch = mpc.branch
    ref, pv, pq = bustypes(bus, gen)

    # generator info
    on_genind = np.where(gen[:,GEN_STATUS])[0]
    gbus = gen[on_genind,GEN_BUS]

    #-----  run the power flow  -----
    t_start = time.time()
    results.iterations = 0
    
    if mpopt.verbose:
        print(' -- DC Power Flow -- \n')
    
    # initial state
    theta0 = bus[:,VA] * (np.pi/180.0)
    
    # build B matrices and phase shift injections
    B, Bf, Pbusinj, Pfinj = makeBdc(bus, branch)
    
    # compute complex bus power injections (generation - load)
    # adjusted for phase shifters and real shunts
    sbus,_ = makeSbus(baseMVA, bus, gen)
    pbus = sbus.real - Pbusinj - bus[:,GS] / baseMVA
    
    # "run" the power flow
    theta,success = dcpf(B, pbus, theta0, ref, pv, pq)
    iterations = 1
    
    # update data matrices with solution
    nl = branch.shape[0]
    branch[:, [QF, QT]] = np.zeros(shape=(nl, 2))
    branch[:, PF] = (Bf * theta + Pfinj) * baseMVA;
    branch[:, PT] = -branch[:, PF]
    bus[:, VM] = np.ones(shape=(bus.shape[0],))
    bus[:, VA] = theta * (180.0/np.pi)
    
    # update Pg for slack generator (1st gen at ref bus)
    # (note: other gens at ref bus are accounted for in Pbus)
    refgen = [i for i,g in enumerate(gbus) if gbus[i] in ref]
    gen[refgen, PG] = gen[refgen, PG] + (B[ref, :] * theta - pbus[ref]) * baseMVA
    
    # -----  output results  -----
    # convert back to original bus numbering & print results
    mpc.bus = bus
    mpc.gen = gen
    mpc.branch = branch
    
    mpc = int2ext(mpc)
    
    # Store the results in a record type
    t_end = time.time()
    results.mpc = mpc
    results.success = success
    results.iterations = iterations
    results.time = t_end-t_start
    
    # zero out result fields of out-of-service gens & branches
    if len(results.mpc.order.gen.status.off)!=0:
        results.mpc.gen[results.mpc.order.gen.status.off, [PG, QG]] = 0
        
    if len(results.mpc.order.branch.status.off)!=0:
        results.mpc.branch[results.mpc.order.branch.status.off, 
                           [PF, QF, PT, QT]] = 0
    if fname != None:
        f = open(os.getcwd()+"/progress/"+fname, 'a')
        text = printpf(results,mpopt)
        f.write(text)
        f.close()
    return results
