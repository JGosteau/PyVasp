#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 10:07:12 2018

@author: gosteau
"""

from lib_gen import *
from lib_poscar import *

s = ['s']
p = ['py', 'pz', 'px']
d = ['dxy' , 'dyz' , 'dz2' , 'dxz' , 'dx2-y2']
f = ['f1' , 'f2' , 'f3', 'f4', 'f5' , 'f6' , 'f7']
ORBS = {'s' : s,'p':p,'d':d,'f':f}

def readfile_to_np(file, skiprows = 0, nrows = None, delimiter = None):
    array = []
    with open(file, 'r') as f :
        for i in range(skiprows):
            f.readline()
        if nrows == None :
            for line in f :
                array.append(line.split(delimiter))
        else : 
            for i in range(nrows) :
                line = f.readline()
                tab = line.split(delimiter)
                if len(tab) != 0 :
                    array.append(line.split(delimiter))
    return np.array(array)

class DOSCAR():    
    def get_def_orbitals(self) :
        if self.lsorbit == True :
            f_orbs = 4
        elif self.spin == 2 :
            f_orbs = 2
        else :
            f_orbs = 1
        n_orbs = int(len(self.occ[0])/f_orbs)
        #print(n_orbs, f_orbs, len(self.occ[0]),self.lorbit)
        if self.lorbit == 11 :
            if n_orbs == 4 :
                orbs = ['s','px','py','pz']
            if n_orbs == 9 :
                orbs = ['s','px','py','pz', 
                        'dxy', 'dyz','dz2', 'dxz','x2-y2']
            if n_orbs == 16 :
                orbs = ['s','px','py','pz', 
                        'dxy', 'dyz','dz2', 'dxz','x2-y2',
                       'f1','f2','f3','f4','f5','f6','f7']
        else :
            if n_orbs == 2 :
                orbs = ['s','p']
            if n_orbs == 3 :
                orbs = ['s','p','d']
            if n_orbs == 4 :
                orbs = ['s','p','d','f']
        return orbs
            
    def get_orbitals(self):
        if self.atom != None :
            orbs = self.get_def_orbitals()
            if self.lsorbit == True :
                spins = ['','_sx','_sy','_sz']
            elif self.spin == 2 :
                spins = ['_up','_dn']
            else :
                spins = ['']
            orbitals = []
            for orb in orbs :
                for s in spins :
                    orbitals.append(orb+s)
            return orbitals
        else :
            if self.spin == 2 :
                orbitals = ['tot_up','tot_dn','sum_up','sum_dn']
            else :
                orbitals = ['tot', 'sum']
            return orbitals

    
    
    def get_ISPIN(self):
        return int(grep("ISPIN", self.path + "OUTCAR")[0][0].split()[2])

    def get_LORBIT(self):
        return int(grep("LORBIT", self.path+"OUTCAR")[0][0].split()[2])

    def get_nDOS(self):
        return int(grep("NEDOS", self.path+"OUTCAR")[0][0].split("=")[1].split()[0])

    def get_LSORBIT(self):
        LSORBIT = grep("LSORBIT", self.path+"OUTCAR")[0][0].split("=")[1].split()[0]
        if LSORBIT == 'T':
            return True
        else :
            return False   
    def get_DOS(self):
        if self.spin == 1 :
            columns = ['Energy', 'dos', 'sum']
        else :
            columns = ['Energy', 'dos_up', 'dos_dn','sum_up', 'sum_dn']
        data = np.float64(readfile_to_np(self.path+"DOSCAR",skiprows = 6, nrows = self.n_dos))
        return data
    
    def get_DOS(self):
        if self.spin == 1 :
            columns = ['Energy', 'dos', 'sum']
        else :
            columns = ['Energy', 'dos_up', 'dos_dn','sum_up', 'sum_dn']
        data = np.float128(readfile_to_np(self.path+"DOSCAR",skiprows = self.ini_dos, nrows = self.n_dos))
        return data
    
    def get_at_DOS(self, atom):
        if self.spin == 1 :
            columns = ['Energy', 'dos', 'sum']
        else :
            columns = ['Energy', 'dos_up', 'dos_dn','sum_up', 'sum_dn']
        tmp = readfile_to_np(self.path+"DOSCAR",skiprows = self.ini_dos+(self.n_dos+1)*atom, nrows = self.n_dos) 
        print(tmp)
        data = np.float128(tmp)
        return data
    
    def __init__(self, path = None, atom = None, ini_dos = 6,ORB_F = False, LORBIT = None):
        
        self.ORB_F = None
        self.atom = []
        self.ini_dos = None
        self.path = None
        self.lorbit = LORBIT
        self.spin = None
        self.n_dos = None
        self.lsorbit = None
        self.energy = 0 
        self.occ = 0
        self.orbitals = None
        if path != None :
            self.ORB_F = ORB_F
            self.atom = [atom]
            self.ini_dos = ini_dos
            self.path = path
            #self.lorbit = self.get_LORBIT()
            self.spin = self.get_ISPIN()
            self.n_dos = self.get_nDOS()
            self.lsorbit = self.get_LSORBIT()
            if atom == None :
                dos = self.get_DOS()
            else : 
                dos = self.get_at_DOS(atom)

            self.energy = dos[:,0]
            self.occ = dos[:,1:]
            self.orbitals = self.get_orbitals()
            
    def __add__(self,other):
        doscar = DOSCAR()
        doscar.atom = []
        for at in self.atom :
            doscar.atom.append(at)
        for at in other.atom :
            doscar.atom.append(at)
        doscar.ini_dos = self.ini_dos
        doscar.path = self.ini_dos
        doscar.lorbit = self.lorbit
        doscar.spin = self.spin
        doscar.n_dos = self.n_dos
        doscar.lsorbit = self.lsorbit
        doscar.energy = self.energy
        doscar.orbitals = self.orbitals
        doscar.energy = self.energy
        doscar.ORB_F = self.ORB_F
        doscar.occ = self.occ + other.occ
        return doscar
    
    def get_orb(self, orb):
        if type(orb) == int :
            return orb
        else :
            index = []
            for l in orb.replace(' ','').split('+') :
                if l not in self.orbitals :
                    try :
                        o,s = l.split('_')
                    except :
                        o = l.split('_')[0]
                        s = ''
                    if o == 'p' :
                        index += self.get_orb('px%s+py%s+pz%s' %(s,s,s))
                    elif o == 'd' :
                        index += self.get_orb('dxy%s+dyz%s+dz2%s+dxz%s+x2-y2%s' %(s,s,s,s,s))
                    elif o == 'f' :
                        index += self.get_orb('f1%s+f2%s+f3%s+f4%s+f5%s+f6%s+f7%s' %(s,s,s,s,s,s,s))
                    elif o == 'tot' :
                        orbs = self.orbitals[0]+s
                        for k in self.orbitals[1:] :
                            orbs+='+'+k+s
                        index += self.get_orb(orbs)
                    else:
                        raise Exception(str(l)+ ' not in ' + str(self.orbitals))
                else :
                    index.append(int(np.where(np.array(self.orbitals) == l)[0][0]))
            return index
    def __getitem__(self, index):
        if type(index) == list or type(index) == np.ndarray :
            res = np.zeros(self.n_dos)
            for i in index :
                res += self[i]
            return res
        else :
            orb = self.get_orb(index)
            if type(orb) == int :
                return self.occ[:,orb]
            else :
                return self[orb]
            
    
     
