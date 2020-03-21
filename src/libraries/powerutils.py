# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 18:48:23 2020

Author: Rounak Meyur
Description: Utilities to run power flow operations for power system analysis
"""
import sys
import time
import numpy as np
from utils import loadcase,bustypes



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
    
    See also
    -------
    rundcpf, runpf
    """
    # arbitrary threshold on |Va| for declaring failure
    theta_threshold = 1e5

    # initialize result vector
    theta = theta0

    # update angles for non-reference buses
    Br_inv = np.linalg.inv(B[pv+pq,:][:,pv+pq])
    pinj = pbus[pv+pq] - (B[pv+pq,ref] * theta0[ref])
    theta([pv+pq]) =  np.dot(Br_inv,pinj)

    # check for presence of *any* warning
    if max(abs(theta)) > theta_threshold:
        print("WARNING!!!: DC power flow failed to solve. Angles above threshold.")
        return theta0
    
    return theta


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


def runpf(casedata, mpopt, fname, solvedcase):
    """
    

    Parameters
    ----------
    casedata : TYPE
        DESCRIPTION.
    mpopt : TYPE
        DESCRIPTION.
    fname : TYPE
        DESCRIPTION.
    solvedcase : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # read data
    mpc = loadcase(casedata)
    
    # get bus index lists of each type of bus
    ref, pv, pq = bustypes(mpc.bus, mpc.gen)

    # generator info
    on_genind = mpc.gen.query('GEN_STATUS > 0').index.tolist()
    gbus = mpc.gen.loc[on_genind,['GEN_BUS']].values

    #-----  run the power flow  -----
    t0 = time.time()
    iterations = 0
    # if mpopt.verbose > 0:
    #     v = mpver('all')
    #     print('\nPYTHONPOWER Version %s, %s', v.Version, v.Date)
    
    if mpopt.verbose > 0:
        print(' -- DC Power Flow -- \n')
    # initial state
    theta0 = mpc.bus['VA'].values * (np.pi/180.0)
    
    # build B matrices and phase shift injections
    [B, Bf, Pbusinj, Pfinj] = makeBdc(baseMVA, bus, branch);
    
    # compute complex bus power injections (generation - load)
    # adjusted for phase shifters and real shunts
    pbus = makeSbus(mpc.baseMVA, mpc.bus, mpc.gen).real - Pbusinj \
        - mpc.bus['GS'].values / mpc.baseMVA
    
    # "run" the power flow
    theta = dcpf(B, pbus, theta0, ref, pv, pq)
    iterations = 1
    
    %% update data matrices with solution
    branch(:, [QF, QT]) = zeros(size(branch, 1), 2);
    branch(:, PF) = (Bf * theta + Pfinj) * baseMVA;
    branch(:, PT) = -branch(:, PF);
    bus(:, VM) = ones(size(bus, 1), 1);
    bus(:, VA) = Va * (180/pi);
    %% update Pg for slack generator (1st gen at ref bus)
    %% (note: other gens at ref bus are accounted for in Pbus)
    %%      Pg = Pinj + Pload + Gs
    %%      newPg = oldPg + newPinj - oldPinj
    refgen = zeros(size(ref));
    for k = 1:length(ref)
        temp = find(gbus == ref(k));
        refgen(k) = on(temp(1));
    end
    gen(refgen, PG) = gen(refgen, PG) + (B(ref, :) * Va - Pbus(ref)) * baseMVA;
    
    
    mpc.et = etime(clock, t0);
    mpc.success = success;
    mpc.iterations = its;
    
    %%-----  output results  -----
    %% convert back to original bus numbering & print results
    [mpc.bus, mpc.gen, mpc.branch] = deal(bus, gen, branch);
    results = int2ext(mpc);
    
    %% zero out result fields of out-of-service gens & branches
    if ~isempty(results.order.gen.status.off)
      results.gen(results.order.gen.status.off, [PG QG]) = 0;
    end
    if ~isempty(results.order.branch.status.off)
      results.branch(results.order.branch.status.off, [PF QF PT QT]) = 0;
    end
    
    if fname
        [fd, msg] = fopen(fname, 'at');
        if fd == -1
            error(msg);
        else
            if mpopt.out.all == 0
                printpf(results, fd, mpoption(mpopt, 'out.all', -1));
            else
                printpf(results, fd, mpopt);
            end
            fclose(fd);
        end
    end
    printpf(results, 1, mpopt);
    
    %% save solved case
    if solvedcase
        savecase(solvedcase, results);
    end
    
    if nargout == 1 || nargout == 2
        MVAbase = results;
        bus = success;
    elseif nargout > 2
        [MVAbase, bus, gen, branch, et] = ...
            deal(results.baseMVA, results.bus, results.gen, results.branch, results.et);
    % else  %% don't define MVAbase, so it doesn't print anything
    end
