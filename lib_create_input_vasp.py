#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 10:50:25 2019

@author: gosteau
"""

import numpy as np


VASP_INCAR_DEF = {'LORBIT' : '.FALSE.'}


HS_ref = {'G' : [0.0,0.0,0.0], 
          'X' : [0.5,0.0,0.0], 'Y' : [0.0,0.5,0.0], 'Z' : [0.0,0.0,0.5],
          'M' : [0.5,0.5,0.0], 'N': [0.0,0.5,0.5], 'O': [0.5,0.0,0.5], 'R' : [0.5,0,0.5],
          'A' : [0.5,0.5,0.5]}
class HS :    
    """
    Classe définissant les points de Haute Symétrie.
    Exemple d'utilisation :
        G = HS('G') # Donne le point HS gamma
        X = HS('X') # Donne le point HS X
        ...
        # la liste des points est défini dans le dictionnaire HS_ref (à completer si on souhaite rajoute des points HS)
        
        # Pour un point de HS symétrie personnalisé de coordonnées kx, ky, kz:
        HS_pers = HS(kx = kx, ky = ky, kz = kz)
        # Il est aussi possible de faire des opérations sur les point de HS
        HS_X_3 = HS('X/3') # Donne aussi le point X/3
        HS_X_3 = HS('X')/3 # Donne aussi le point X/3
        Le dictionnaire HS_points donne la liste des principaux points de HS.
        G = HS_points('G')
        X = HS_points('X')
        X_3 = HS_points('X')/3
    """
    def __init__(self,name= None, kx=0.0, ky = 0.0,kz = 0.0):
        try :
            sp = name.split('/')
            hs = HS(sp[0])/float(sp[1])
            self.__dict__ = hs.__dict__
        except :
            if name in HS_ref.keys():
                self.name = name
                self.kx = HS_ref[name][0]
                self.ky = HS_ref[name][1]
                self.kz = HS_ref[name][2]
            else :
                self.name = name
                self.kx = kx
                self.ky = ky
                self.kz = kz
    
    def __div__(self,div):
        name = str(self.name)+'/'+str(div)
        return HS(name, self.kx/div, self.ky/div, self.kz/div)   

    def __mul__(self,mul):
        name = str(self.name)+'*'+str(mul)
        return HS(name, self.kx*mul, self.ky*mul, self.kz*mul)        
    
    def __repr__(self):
        return "%s : %f %f %f" %(self.name, self.kx, self.ky, self.kz)
    
    def __getitem__(self, index):
        if type(index) == int:
            index = get_k(index)
        if index == "kx" :
            return self.kx
        elif index == "ky" :
            return self.ky
        elif index == "kz" :
            return self.kz
        
    def __setitem__(self, index, val):
        if type(index) == int:
            index = get_k(index)
        if index == "kx" :
            self.kx = val
        elif index == "ky" :
            self.ky = val
        elif index == "kz" :
            self.kz = val
    
    def __add__(self, other):
        hs = HS()
        hs.name = self.name+'+'+other.name
        for i in ['kx','ky','kz']:
            hs[i] = self[i] + other[i]
        return hs
    
    def __sub__(self, other):
        hs = HS()
        hs.name = self.name+'-'+other.name
        for i in ['kx','ky','kz']:
            hs[i] = self[i] - other[i]
        return hs
    
    def __eq__(self,other):
        if self['kx'] == other['kx'] and self['ky'] == other['ky'] and self['kz'] == other['kz'] :
            return True
        else :
            return False
        

G = HS('G')
X = HS('X')
M = HS('M')
Z = HS('Z')
A = HS('A')
R = HS('R')
HS_points = {'G' : G, 'X' : X, 'M' : M, 'Z' : Z, 'A' : A, 'R' : R}

def_HS = [[G,X],
          [X,M],
          [M,G],
          [G,Z],
          [Z,R],
          [R,A],
          [A,Z]]



class KPOINTS:
    def __init__(self, Type, nkpt,line):
        self.type = Type
        self.nkpt = nkpt
        self.line = line
        
    def __repr__(self):
        string = ''
        if self.type[0] == 'G' or self.type[0] == 'M':
            if type(self.nkpt) == np.ndarray or type(self.nkpt) == list :
                nx, ny, nz = self.nkpt
            else :
                nx = self.nkpt; ny = self.nkpt; nz = self.nkpt
            string += '%dx%dx%d\n0\n' % (nx,ny,nz)
            string += '%s\n' %(self.type)
            string += '%d %d %d\n' %(nx,ny,nz)
            string += '0 0 0'
            
        elif self.type == 'line':
            
            string += 'BS calculation\n %d\n' % (self.nkpt)
            #string += '%s\n' %(self.type)
            string += 'Line-mode\n'
            string += 'rec\n'
            for kpt in np.linspace(self.line[0],self.line[1],self.nkpt):
                string += " %.8f   %.8f   %.8f\n" %(kpt, self.line[0],0.5)
                string += " %.8f   %.8f   %.8f\n\n" %(kpt, self.line[1],0.5)
                
            
        elif self.type == 'From-line':
            for hs_line in self.line :
                BOOL = False
                hs1 = hs_line[0]
                hs2 = hs_line[1]
                hs = hs2-hs1
                tab = []
                for i in np.arange(self.nkpt):
                    hs_i = hs1+hs*(i/(self.nkpt-1))
                    string += (" %.8f   %.8f   %.8f   0\n"%(hs_i['kx'],hs_i['ky'],hs_i['kz']))
                    tab.append(hs_i)
                
            
        elif self.type == 'FS':
            string += 'FS calculation\n %d\n' % (self.nkpt)
            #string += '%s\n' %(self.type)
            string += 'Line-mode\n'
            string += 'rec\n'
            for i in np.arange(self.nkpt):
                point_1 = [-0.5,i/(self.nkpt-1)-0.5,0.5]
                point_2 = [0.5,i/(self.nkpt-1)-0.5,0.5]
                string += (" %.8f   %.8f   %.8f\n"%(point_1[0], point_1[1], point_1[2]))
                string += (" %.8f   %.8f   %.8f\n\n"%(point_2[0], point_2[1], point_2[2]))
        return string


class INCAR:
    def add(self, comment = False, **kwargs):
        for key in kwargs.keys():
            if comment :
                self.dict['#'+key] = kwargs.get(key)
            else :
                self.dict[key] = kwargs.get(key)
            
    def __init__(self,**kwargs):
        self.dict = {}
        self.add(**kwargs)
        
    def __getitem__(self, key):
        if key in self.dict.keys() :
            return self.dict[key]
        elif key in VASP_INCAR_DEF.keys() :
            string = "%s not defined in this incar. Vasp will take the default value :\n%s = %s" %(key, key, VASP_INCAR_DEF[key])
            print(string)
            return None
        else :
            string = "%s parameter does not exist" %(key)
            print(string)
            return None
    
    
    def __repr__(self):
        string = ""
        for key in self.dict.keys():
            var = self.dict[key]
            if type(var) == type(None) :
                string += "%s\n" %(key)
            else  :
                if type(var) == bool:
                    if var == True :
                        variable = '.TRUE.'
                    elif var == False :
                        variable = '.FALSE.'
                elif type(var) == list or type(var) == np.ndarray :
                    variable = ''
                    for i in var :
                        variable += str(i) + ' '
                else :
                    variable = str(var)
                string += "%s = %s\n" %(key, variable)
        return string