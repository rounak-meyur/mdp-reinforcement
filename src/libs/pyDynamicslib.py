# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:10:32 2019

@author: rounak
"""

import numpy as np
import scipy.sparse as sp

#%% Global Variables
global base_mva, base_kV, omega, Zbase
base_mva = 100.0
base_kV = 33.0
omega = 2*np.pi*60
Zbase = (base_kV**2)/base_mva


#%% Functions
def KronRed(A,col1,col2):
    '''
    '''
    K1 = A[col1,col1.T]
    K2 = A[col1,col2.T]
    K3 = A[col2,col1.T]
    K4 = A[col2,col2.T]
    return K1 - np.matmul(np.matmul(K2,np.linalg.inv(K4)),K3)


def SlowCoherency(A,n):
    '''
    Identifies the largest n eigen values of the system matrix A and returns
    them along with the corresponding eigen vectors.
    
    Input:  A: System matrix
            n: number of areas (slow coherent areas)
    Output: w: largest n eigen values
            v: corresponding eigenvectors
    '''
    w,v = np.linalg.eig(A)
    idx = w.argsort()[::-1]
    w = w[idx[0:n]]
    v = v[:,idx[0:n]]
    return w,v


def Group(V):
    '''
    '''
    nmach,narea = np.shape(V)
    order = np.array(range(nmach))
    Vaug = np.row_stack((np.zeros(shape=(1,narea+1)),
                         np.column_stack((np.zeros(shape=(nmach,1)),V))))
    # Order the generators to determine the reference buses 
    for i in range(narea):
        Vaug = Vaug[1:nmach-i+1,1:narea-i+1]
        valmax = np.max(abs(Vaug),axis=0)
        rowind = np.argmax(abs(Vaug),axis=0)
        colref = np.argmax(valmax)
        rowref = rowind[colref]
        # Permutations
        order[[i,rowref]] = order[[rowref,i]]
        Vaug[[0,rowref],:] = Vaug[[rowref,0],:]
        Vaug[:,[0,colref]] = Vaug[:,[colref,0]]
        # Gaussian Elimination
        v = (Vaug[:,0]/Vaug[0,0]).reshape(-1,1)
        X = np.matmul(v,Vaug[0,:].reshape(1,-1))
        Vaug = Vaug - X
    
    # Build the modified eigenvector matrix and get area assignment
    Vnew = V[order,:]
    V1 = Vnew[:narea,:]
    V2 = Vnew[narea:,:]
    L = np.matmul(V2,np.linalg.inv(V1))
    k = np.concatenate((np.array(range(narea)),np.argmax(L,axis=1)))
    return order+1,k


class SmallSignal:
    """
    """
    def __init__(self,bus,line,tx,mach,pfsol):
        '''
        '''
        gen = np.array(range(len(bus),len(bus)+len(mach))).reshape(-1,1)
        nongen = np.array(range(len(bus))).reshape(-1,1)
        YbusDC = self.__GetYbusDC(bus,line,tx,mach)
        self.Ybusred = KronRed(YbusDC,gen,nongen)
        (self.E,self.delta) = self.__IntVoltage(bus,mach,pfsol)
        self.Asys = np.zeros(shape=(2*len(mach),2*len(mach)),dtype=float)
        return
        
        
    def __IntVoltage(self,bus,mach,pfsol):
        '''
        Evaluates the generator bus internal voltage magnitudes and angles.
        '''
        gbus = np.array([bus.ind(m['number']) for m in mach]).reshape(-1,1)
        cos_theta = np.cos(pfsol.theta[gbus,0])
        sin_theta = np.sin(pfsol.theta[gbus,0])
        Sgen = pfsol.Pgen[gbus,0] + 1j*pfsol.Qgen[gbus,0]
        V = pfsol.V[gbus,0]*(cos_theta+1j*sin_theta)
        I = (Sgen/V).conjugate()
        X = np.array([m['xd'] for m in mach]).reshape(-1,1)*(1j)
        E = V + (I*X)
        
        Emag = abs(E); Eang = np.arctan2(E.imag,E.real)
        return (Emag,Eang)
    
    
    def __GetYbusDC(self,bus,line,tx,mach):
        '''
        Evaluates the Ybus matrix under DC approximation for a given power 
        system network. Ignores the resistances and shunt admittances.
        '''
        Ybus = np.zeros(shape=(len(bus)+len(mach),len(bus)+len(mach)),\
                        dtype=complex)
        # Line reactances
        for k in range(len(line)):
            # find index of buses
            ind_fbus = bus.ind(line[k]['frombus'])
            ind_tbus = bus.ind(line[k]['tobus'])
            # evaluate reactance
            x = line[k]['x']*line[k]['length']*base_mva/line[k]['mva']
            # add off-diagonal elements of Ybus
            Ybus[ind_fbus,ind_tbus] += -1/x
            Ybus[ind_tbus,ind_fbus] += -1/x
        # Transformer reactances
        for k in range(len(tx)):
            # find index of buses
            ind_fbus = bus.ind(tx[k]['frombus'])
            ind_tbus = bus.ind(tx[k]['tobus'])
            # evaluate reactance
            x = tx[k]['x']*base_mva/tx[k]['mva']
            # add off-diagonal elements of Ybus
            Ybus[ind_fbus,ind_tbus] += -1/x
            Ybus[ind_tbus,ind_fbus] += -1/x
        for m in range(len(mach)):
            # find index of machine bus
            ind_mbus = bus.ind(mach[m]['number'])
            x = mach[m]['xd']*base_mva/mach[m]['rate']
            Ybus[ind_mbus,len(bus)+m] = -1/x
            Ybus[len(bus)+m,ind_mbus] = -1/x
        
        for n in range(len(bus)+len(mach)):
            Ybus[n,n] = -np.sum(Ybus[n,:])
        
        return Ybus
    
    
    def BuildStateMatrix(self,mach):
        '''
        '''
        ngen = len(mach)
        Ks = (self.E.dot(self.E.T))*np.cos(self.delta-self.delta.T)*self.Ybusred
        #Ks = self.Ybusred
        Kd = np.identity(ngen)*10
        Minv = np.diag([1/((2*m['h'])*(m['rate']/base_mva)) for m in mach])
        
        A1 = np.zeros(shape=(ngen,ngen))
        A2 = omega*np.identity(ngen)
        A3 = -Minv.dot(Ks)
        A4 = -Minv.dot(Kd)
        
        self.Asys = np.concatenate((np.concatenate((A1,A2),axis=1),\
                                    np.concatenate((A3,A4),axis=1)),axis=0)
        return
    
    

class Machine:
    """
    Generates dynamic models for a electromechanical machine model with the 
    following paramters.
    """
    def __init__(self,mach,bus):
        '''
        '''
        # Machine data/parameters
        self.__ngen = len(mach)
        self.busind = [mach.ind(mach[i]['number'],mach[i]['id']) \
                       for i in range(self.__ngen)]
        self.scale = np.array([base_mva/mach[i]['rate'] \
                                    for i in range(len(mach))]).reshape(-1,1)
        self.xd = np.array([mach[i]['xd'] \
                            for i in range(self.__ngen)]).reshape(-1,1)
        self.D = np.array([mach[i]['d'] \
                           for i in range(self.__ngen)]).reshape(-1,1)
        self.H = np.array([mach[i]['h'] \
                           for i in range(self.__ngen)]).reshape(-1,1)
        
        # Get data from solved power flow
        Vmag = np.array([bus[i]['pu'] for i in self.busind]).reshape(-1,1)
        theta = np.array([bus[i]['angle'] \
                          for i in self.busind]).reshape(-1,1)*np.pi/180
        V = Vmag*np.exp(1j*(theta))
        # Real/Reactive power on system MVA base
        self.P = np.array([bus[i]['pgen'] \
                               for i in self.busind]).reshape(-1,1)/base_mva
        self.Q = np.array([bus[i]['qgen'] \
                               for i in self.busind]).reshape(-1,1)/base_mva
        
        # Mechanical power on generator mva base
        self.Pmech = self.P*self.scale
        
        # Generator current on generator MVA base
        I = self.scale*((self.P+1j*self.Q).conjugate()/V)
       
        # Computing stator variables
        self.psi = V + 1j*self.xd*I
        
        # Computing stator variables in rotating frame of reference
        self.delta = np.angle(self.psi)
        self.speed = np.ones(shape=(self.__ngen,1))
        rot = 1j*np.exp(-1j*self.delta)
        Erot = self.psi*rot; Irot = I*rot; Vrot = V*rot
        self.Ed = Erot.real; self.Eq = Erot.imag
        self.Vd = Vrot.real; self.Vq = Vrot.imag
        self.Idg = Irot.real; self.Iqg = Irot.imag
        self.Vex = self.Eq
        
        # Generator currents on system mva base
        self.Id = self.Idg/self.scale
        self.Iq = self.Iqg/self.scale
        
        # Initialize delta variables
        self.d_delta = np.zeros(shape=(self.__ngen,1))
        self.d_speed = np.zeros(shape=(self.__ngen,1))
        self.A = np.zeros(shape=(2*self.__ngen,2*self.__ngen))
        return
    
    
    def GenDyn(self,ref,Yred):
        '''
        '''
        self.d_delta = np.zeros(shape=(self.__ngen,1))
        self.d_speed = np.zeros(shape=(self.__ngen,1))
        
        self.delta = self.delta - ref
        self.psi = (np.sin(self.delta)*self.Ed \
                    + np.cos(self.delta)*self.Eq) \
                    + 1j*(-np.cos(self.delta)*self.Ed \
                          + np.sin(self.delta)*self.Eq)
        
        # Calculation on system base
        I = Yred.dot(self.psi)
        self.Id = I.real*np.sin(self.delta) - I.imag*np.cos(self.delta)
        self.Iq = I.real*np.cos(self.delta) + I.imag*np.sin(self.delta)
        
        # Calculation on generator mva base
        self.Idg = self.Id*self.scale; self.Iqg = self.Iq*self.scale
        self.Vd = self.Ed + (self.xd*self.Iqg)
        self.Vq = self.Eq - (self.xd*self.Idg)
        
        # Calculation of real/reactive power on system base
        self.P = self.Eq*self.Iq + self.Ed*self.Id
        self.Q = self.Eq*self.Id - self.Ed*self.Iq
        
        self.d_delta = omega*(self.speed-np.ones(shape=(self.__ngen,1)))
        self.d_speed = (self.Pmech-(self.P*self.scale)-\
                        (self.D*(self.speed-np.ones(shape=(self.__ngen,1)))))\
                        /(2*self.H)
        return
    
    
    def Perturb(self,Yred,perturb = 1e-3,ref=0):
        '''
        '''
        for i in range(self.__ngen):
            pert = np.zeros(shape=(self.__ngen,1))
            nominal = self.delta
            pert[i,0] = perturb * self.delta[i,0]
            self.delta = nominal + pert
            self.GenDyn(ref,Yred)
            self.A[:self.__ngen,i] = self.d_delta[:,0]/pert[i,0]
            self.A[self.__ngen:,i] = self.d_speed[:,0]/pert[i,0]
            self.delta = nominal
        for i in range(self.__ngen):
            pert = np.zeros(shape=(self.__ngen,1))
            nominal = self.speed
            pert[i,0] = perturb * self.speed[i,0]
            self.speed = nominal + pert
            self.GenDyn(ref,Yred)
            self.A[:self.__ngen,i+self.__ngen] = self.d_delta[:,0]/pert[i,0]
            self.A[self.__ngen:,i+self.__ngen] = self.d_speed[:,0]/pert[i,0]
            self.speed = nominal
        return self.A[self.__ngen:,:self.__ngen]






















