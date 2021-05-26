# -*- coding: utf-8 -*-
"""
Created on Mon Dec 14 23:09:13 2020

@author: rouna
"""
import numpy as np
global base_mva, base_kV, omega, Zbase
base_mva = 100.0
base_kV = 33.0
omega = 2*np.pi*60
Zbase = (base_kV**2)/base_mva

#%% Classes
class iBus(list):
    """
    A class with attributes and methods for a single bus instance.
    The attributes corresponding to each bus instance are as follows:
        0: bus number (integer type)
        1: bus type (integer type, 1: swing bus, 2: generator bus, 3: load bus)
        2: per unit voltage magnitude at the bus (floating type)
        3: bus voltage angle in degrees (floating type)
        4: real power generated or injected at the bus (floating type)
        5: reactive power generated or injected at the bus (floating type)
        6: real power load or consumption at the bus (floating type)
        7: reactive power load or consumption at the bus (floating type)
        8: Minimum reactive power generation at the bus (floating type)
        9: Maximum reactive power generation at the bus (floating type)
        10: Shunt conductance at the bus in MW (floating type)
        11: Shunt susceptance at the bus in MVAR (floating type)
        12: Base voltage of the bus in pu (floating type)
    The attribute values can be accessed using 'list' methods. For example, if
    b is an instance of the class iBus, b[0]=b['number'] will access the
    bus number of the bus.
    """
    def __init__(self,aList):
        '''
        Constructor to initiate a bus instance/object.
        A bus object has two kinds of data. The first data refers a number to
        each variable name (for example bus number, type etc.) and the second 
        data is the input data in the form of list. These two sets of data
        would be used to map the input values to corresponding variable name.
        '''
        super(iBus,self).__init__(aList)
        self.__dict__['d'] = {'number':0, 'type':1, 'pu':2, 'angle':3, 'pgen':4,
                     'qgen':5, 'pload':6, 'qload':7, 'qmin':8, 'qmax':9, 'g':10,
                     'b':11, 'kv':12}
        self.__dict__['data'] = aList
    
    def __getitem__(self,key):
        '''
        A method to get item from class using the key defined in the constructor
        '''
        if isinstance(key,str):
            return self.data[self.d[key.lower()]]
        else:
            return self.data[key]
    
    def __setitem__(self,key,value):
        '''
        A method to get set value for a particular attribute of the instance.
        '''
        if isinstance(key,str):
            self.data[self.d[key.lower()]] = value
        else:
            self.data[key] = value


class iLine(list):
    """
    A class with attributes and methods for a single line instance. The 
    attributes corresponding to each line instance are as follows:
        0: from bus number (integer type)
        1: to bus number (integer type)
        2: line ID (string type)
        3: resistance of line in pu/km (floating type)
        4: reactance of line in pu/km (floating type)
        5: shunt susceptance of line in pu/km (floating type)
        6: status of line (binary type)
        7: length of line in km (floating type)
        8: MVA rating of line (floating type)
        9: kV rating of the line (floating type)
    The attribute values can be accessed using 'list' methods. For example, if
    b is an instance of the class iLine, b[2]=b['id'] will access the
    line ID of the line.
    """
    def __init__(self,aList):
        '''
        Constructor to initiate a line instance/object.
        A line object has two kinds of data. The first data refers a number to
        each variable name (for example from bus, to bus etc.) and the second 
        data is the input data in the form of list. These two sets of data
        would be used to map the input values to corresponding variable name.
        '''
        super(iLine,self).__init__(aList)
        self.__dict__['d'] = {'frombus':0,'tobus':1,'id':2,'r':3,'x':4,'b':5,
                     'st':6, 'length':7, 'mva':8, 'kv':9}
        self.__dict__['data'] = aList
    
    def __getitem__(self,key):
        '''
        A method to get item from class using the key defined in the constructor
        '''
        if isinstance(key,str):
            return self.data[self.d[key.lower()]]
        else:
            return self.data[key]
    
    def __setitem__(self,key,value):
        '''
        A method to get set value for a particular attribute of the instance.
        '''
        if isinstance(key,str):
            self.data[self.d[key.lower()]] = value
        else:
            self.data[key] = value


class iTx(list):
    """
    A class with attributes and methods for a single two winding transformer
    instance. The attributes corresponding to each transformer instance are as 
    follows:
        0: from bus number (integer type)
        1: to bus number (integer type)
        2: transformer ID (string type)
        3: resistance of transformer in pu (floating type)
        4: reactance of transformer in pu (floating type)
        5: shunt conductance of transformer in pu (floating type)
        6: shunt susceptance of transformer in pu (floating type)
        7: tap ratio of transformer (floating type)
        8: status of transformer (binary type)
        9: MVA rating of transformer (floating type)
        10: kV rating of HV side of transformer (floating type)
        11: kV rating of LV side of transformer (floating type)
    The attribute values can be accessed using 'list' methods. For example, if
    b is an instance of the class iTx, b[2]=b['id'] will access the
    transformer ID of the transformer.
    """
    def __init__(self,aList):
        '''
        Constructor to initiate a transformer instance/object.
        A transformer object has two kinds of data. The first data refers a number 
        to each variable name (for example from bus, to bus etc.) and the second 
        data is the input data in the form of list. These two sets of data
        would be used to map the input values to corresponding variable name.
        '''
        super(iTx,self).__init__(aList)
        self.__dict__['d'] = {'frombus':0,'tobus':1,'id':2,'r':3,'x':4,'g':5,
                     'b':6, 'tap':7, 'status':8, 'mva':9, 'hv':10, 'lv':11}
        self.__dict__['data'] = aList
    
    def __getitem__(self,key):
        '''
        A method to get item from class using the key defined in the constructor
        '''
        if isinstance(key,str):
            return self.data[self.d[key.lower()]]
        else:
            return self.data[key]
    
    def __setitem__(self,key,value):
        '''
        A method to get set value for a particular attribute of the instance.
        '''
        if isinstance(key,str):
            self.data[self.d[key.lower()]] = value
        else:
            self.data[key] = value


class iMach(list):
    """
    A class with attributes and methods for a single machine instance.
    The attributes corresponding to each machine instance are as follows:
        0: generator/machine ID (string type)
        1: bus number (integer type)
        2: rating of machine (floating type)
        3: sub-transient reactance (floating type)
        4: machine inertia (floating type)
    The attribute values can be accessed using 'list' methods. For example, if
    g is an instance of the class iMach, g[4]=b['H'] will access the machine
    inertia of the generator.
    """
    def __init__(self,aList):
        '''
        Constructor to initiate a machine or generator instance/object.
        A machine object has two kinds of data. The first data refers a number
        to each variable name (for example bus number, inertia etc.) and the
        second data is the input data in the form of list. These two sets of data
        would be used to map the input values to corresponding variable name.
        '''
        super(iMach,self).__init__(aList)
        self.__dict__['d'] = {'number':0, 'id':1, 'rate':2, 'xd':3, 'h':4,
                     'd':5}
        self.__dict__['data'] = aList
    
    def __getitem__(self,key):
        '''
        A method to get item from class using the key defined in the constructor
        '''
        if isinstance(key,str):
            return self.data[self.d[key.lower()]]
        else:
            return self.data[key]
    
    def __setitem__(self,key,value):
        '''
        A method to get set value for a particular attribute of the instance.
        '''
        if isinstance(key,str):
            self.data[self.d[key.lower()]] = value
        else:
            self.data[key] = value


class Bus(list):
    """
    A class to instantiate all the buses in the system of interest
    """
    def __getitem__(self,key):
        '''
        Returns the bus instance corresponding to key
        '''
        return iBus(self.data[key])
    
    def __iter__(self):
        '''
        Returns as class objects during iterations
        '''
        for p in self.data:
            yield iBus(p)
        
    def __init__(self,*args):
        '''
        '''
        if len(args)==1 and isinstance(args[0],str):
            # Initialize through a csv file
            f = open(args[0],'r')
            cdata = [line.strip('\n').split(',') for line in f.readlines()]
            f.close()
            # Append each line of data one at a time
            self.data = []
            dicref = []
            for k in range(1,len(cdata)):
                L = [int(cdata[k][0]),int(cdata[k][1]),float(cdata[k][2]),
                     float(cdata[k][3]),float(cdata[k][4]),float(cdata[k][5]),
                     float(cdata[k][6]),float(cdata[k][7]),float(cdata[k][8]),
                     float(cdata[k][9]),float(cdata[k][10]),float(cdata[k][11]),
                     float(cdata[k][12])]
                self.__dict__['data'].append(L)
                dicref.append(L[0])
            super(Bus,self).__init__(self.__dict__['data'])
            self.identify = dict(zip(dicref,range(len(cdata))))
        elif len(args) == 1:
            # Initialize through a single list of data
            self.__dict__['data'] = args[0]
            super(Bus,self).__init__(self.__dict__['data'])
            dicref = [col[0] for col in args[0]]
            self.identify = dict(zip(dicref,range(len(args[0]))))
        else:
            # Exception Handling
            print("Wrong number/type of arguments:",str(len(args)))
            print("To get data from csv file use: '%pathname%filename.csv'")
            print('To instantiate use a list as the only argument')
    
    def ind(self,number):
        '''
        Identifies the index of the bus in the bus list from the bus number 
        provided as input.
        '''
        ret = []
        key = int(number)
        if key in self.identify: ret = self.identify[key]
        return ret
    
    def update(self,pfsol):
        '''
        Updates the bus data after a power flow operation.
        '''
        for k in range(len(self.data)):
            if self.data[k][1] == 1:
                self.data[k][4] = pfsol.Pgen[k,0]*base_mva
                self.data[k][5] = pfsol.Qgen[k,0]*base_mva
            elif self.data[k][1] == 2:
                self.data[k][2] = pfsol.V[k,0]
                self.data[k][3] = pfsol.theta[k,0]*(180/np.pi)
                self.data[k][5] = pfsol.Qgen[k,0]*base_mva
            elif self.data[k][1] == 3:
                self.data[k][2] = pfsol.V[k,0]
                self.data[k][3] = pfsol.theta[k,0]*(180/np.pi)
        return



class Line(list):
    """
    A class to instantiate all the lines in the system of interest.
    """
    def __getitem__(self,key):
        '''
        Returns the line instance corresponding to key
        '''
        return iLine(self.data[key])
        
    def __iter__(self):
        '''
        Returns as class objects during iterations
        '''
        for p in self.data:
            yield iLine(p)
            
    def __init__(self,*args):
        '''
        '''
        if len(args)==1 and isinstance(args[0],str):
            # Initialize through a csv file
            f = open(args[0],'r')
            cdata = [line.strip('\n').split(',') for line in f.readlines()]
            f.close()
            # Append each line of data one at a time
            self.data = []
            dicref = []
            for k in range(1,len(cdata)):
                L = [int(cdata[k][0]),int(cdata[k][1]),str(cdata[k][2]).strip(),
                     float(cdata[k][3]),float(cdata[k][4]),float(cdata[k][5]),
                     int(cdata[k][6]),float(cdata[k][7]),float(cdata[k][8]),
                     float(cdata[k][9])]
                self.__dict__['data'].append(L)
                dicref.append((L[0],L[1],L[2]))
            super(Line,self).__init__(self.__dict__['data'])
            self.identify = dict(zip(dicref,range(len(cdata))))
        elif len(args) == 1:
            # Initialize through a single list of data
            self.__dict__['data'] = args[0]
            super(Line,self).__init__(self.__dict__['data'])
            dicref = [(col[0],col[1],col[2].strip()) for col in args[0]]
            self.identify = dict(zip(dicref,range(len(args[0]))))
        else:
            # Exception Handling
            print("Wrong number/type of arguments:",str(len(args)))
            print("To get data from csv file use: '%pathname%filename.csv'")
            print('To instantiate use a list as the only argument')
    
    def ind(self,fbus,tbus,bid):
        '''
        Identifies the index of the line in the line list from the from bus,
        to bus and line id provided as input.
        '''
        ret = []
        key = (int(fbus),int(tbus),str(bid).strip())
        if key in self.identify: ret = self.identify[key]
        return ret
            

class Tx(list):
    """
    A class to instantiate all the two winding transformers in the system of 
    interest.
    """
    def __getitem__(self,key):
        '''
        Returns the transformer instance corresponding to key
        '''
        return iTx(self.data[key])
        
    def __iter__(self):
        '''
        Returns as class objects during iterations
        '''
        for p in self.data:
            yield iTx(p)
            
    def __init__(self,*args):
        '''
        '''
        if len(args)==1 and isinstance(args[0],str):
            # Initialize through a csv file
            f = open(args[0],'r')
            cdata = [line.strip('\n').split(',') for line in f.readlines()]
            f.close()
            # Append each line of data one at a time
            self.data = []
            dicref = []
            for k in range(1,len(cdata)):
                L = [int(cdata[k][0]),int(cdata[k][1]),str(cdata[k][2]).strip(),
                     float(cdata[k][3]),float(cdata[k][4]),float(cdata[k][5]),
                     float(cdata[k][6]),float(cdata[k][7]),int(cdata[k][8]),
                     float(cdata[k][9]),float(cdata[k][10]),float(cdata[k][11])]
                self.__dict__['data'].append(L)
                dicref.append((L[0],L[1],L[2]))
            super(Tx,self).__init__(self.__dict__['data'])
            self.identify = dict(zip(dicref,range(len(cdata))))
        elif len(args) == 1:
            # Initialize through a single list of data
            self.__dict__['data'] = args[0]
            super(Tx,self).__init__(self.__dict__['data'])
            dicref = [(col[0],col[1],col[2].strip()) for col in args[0]]
            self.identify = dict(zip(dicref,range(len(args[0]))))
        else:
            # Exception Handling
            print("Wrong number/type of arguments:",str(len(args)))
            print("To get data from csv file use: '%pathname%filename.csv'")
            print('To instantiate use a list as the only argument')
    
    def ind(self,fbus,tbus,bid):
        '''
        Identifies the index of the transformer in the transformer list from 
        the from bus, to bus and transformer id provided as input.
        '''
        ret = []
        key = (int(fbus),int(tbus),str(bid).strip())
        if key in self.identify: ret = self.identify[key]
        return ret


class Mach(list):
    """
    A class to instantiate all the machines in the system of interest
    """
    def __getitem__(self,key):
        '''
        Returns the bus instance corresponding to key
        '''
        return iMach(self.data[key])
    
    def __iter__(self):
        '''
        Returns as class objects during iterations
        '''
        for p in self.data:
            yield iMach(p)
            
    def __init__(self,*args):
        '''
        '''
        if len(args)==1 and isinstance(args[0],str):
            # Initialize through a csv file
            f = open(args[0],'r')
            cdata = [line.strip('\n').split(',') for line in f.readlines()]
            f.close()
            # Append each line of data one at a time
            self.data = []
            dicref = []
            for k in range(1,len(cdata)):
                L = [int(cdata[k][1]),str(cdata[k][0]),float(cdata[k][2]),
                     float(cdata[k][3]),float(cdata[k][4]),float(cdata[k][5])]
                self.__dict__['data'].append(L)
                dicref.append((L[0],L[1]))
            super(Mach,self).__init__(self.__dict__['data'])
            self.identify = dict(zip(dicref,range(len(cdata))))
        elif len(args) == 1:
            # Initialize through a single list of data
            self.__dict__['data'] = args[0]
            super(Mach,self).__init__(self.__dict__['data'])
            dicref = [(col[0],col[1]) for col in args[0]]
            self.identify = dict(zip(dicref,range(len(args[0]))))
        else:
            # Exception Handling
            print("Wrong number/type of arguments:",str(len(args)))
            print("To get data from csv file use: '%pathname%filename.csv'")
            print('To instantiate use a list as the only argument')

    def ind(self,number,gid):
        '''
        Identifies the index of the machine in the machine list from the bus 
        number and generator id provided as input.
        '''
        ret = []
        key = (int(number),str(gid))
        if key in self.identify: ret = self.identify[key]
        return ret