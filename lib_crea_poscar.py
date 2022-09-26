#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:56:41 2019

@author: gosteau
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

AXIS = {'x' : 0, 'y' : 1, 'z' : 2}

class POSCAR():
    def __init__(self):
        self.name = 'POSCAR'
        self.scale = 1
        self.axis  = np.array([[1.,0.,0.],
                               [0.,1.,0.],
                               [0.,0.,1.]])
        self.pos = pd.DataFrame(columns = ['at','x','y','z'])
        self.atoms = {}
        self.repr = "Direct"
        
    def add_atom(self, at,x,y,z):
        self.pos = pd.DataFrame.append(self.pos,{'at': at, 'x' : x, 'y' : y,'z' : z}, ignore_index = True)
        if at in self.atoms.keys():
            self.atoms[at] += 1
        else :
            self.atoms[at] = 1
      
    def add_atoms(self,atoms):
        for i in atoms:
            at,x,y,z = i
            self.add_atom(at,x,y,z)
            
    def set_lattice(self, Vacuum = [False,False,False]):
        for ax in AXIS :
            i = AXIS[ax]
            v = Vacuum[i]
            if v != False :
                self.axis[i][i] = self.pos[ax].max()-self.pos[ax].min() + v
            else :
                self.axis[i][i] = self.pos[ax].max()
    
    def __repr__(self):
        string = ''
        string += "%s\n" %(self.name)
        string += "   %.14f\n" %(self.scale)
        
        for i in self.axis:
            string += ' '
            for j in i :
                string += '    %.14f' %(j)
            string += '\n'
        string += "   "
        string2 = "   "
        for at in self.atoms.keys():
            string += "%3s   " %(at)
            string2 += "%3d   "%(self.atoms[at])
        string += "\n" + string2+ "\n"
        string += "%s\n" %(self.repr)
        
        for at in self.atoms.keys() :
            for i in self.pos.loc[self.pos['at'] == at].index:
                for ax in AXIS.keys():
                    j = AXIS[ax]
                    norm = np.sqrt(self.axis[j][0]**2+self.axis[j][1]**2+self.axis[j][2]**2)
                    string += "   %.14f" %(self.pos[ax][i]/norm)
                string +='\n'
        return string
        
        
        
poscar = POSCAR()
"""
poscar.add_atom('Sr', 0,0,0)
poscar.add_atom('Ti', 0.5,0.5,0.5)
poscar.add_atom('O', 0.5,0.5,0)
poscar.add_atom('O', 0.5,0,0.5)
poscar.add_atom('O', 0,0.5,0.5)
poscar.add_atom('Sr', 0,0,1)
"""     
ax = 3.905
ay = ax
az = ax
v = 15
num = 7
"""
poscar.name = '%.1f SrTiO3' %(num)
poscar.atoms = {'Sr':0,'Ti':0,'O':0}

for i in range(int(num*2)):
    z = i/2*az+v/2
    if i%2 == 1:
        ATOMS = [['Sr',0,0,z],
                 ['O',ax/2,ay/2,z]]
    else :
        ATOMS = [['Ti',ax/2,ay/2,z],
                 ['O',0,ay/2,z],
                 ['O',ax/2,0,z]]
    poscar.add_atoms(ATOMS)

poscar.set_lattice([False,False,v])
poscar.axis[0]*=2
poscar.axis[1]*=2
"""
poscar.atoms = {'Sr':0,'Ti':0,'O':0}

ax = 3.905
ay = ax
az = ax
v = 15
num = 7

poscar.name = '%.1f STO' %(num)
for i in range(int(num)):
    z = i*az+v/2
    ATOMS = [['Sr',0,0,z+az/2],
             ['O',ax/2,ay/2,z+az/2],
             ['Ti',ax/2,ay/2,z],
             ['O',ax/2,0,z],
             ['O',0,ay/2,z]]
    poscar.add_atoms(ATOMS)

poscar.set_lattice([False,False,v])
poscar.axis[0][0]=ax
poscar.axis[1][1]=ay
poscar.axis[2][2]=az*num+v
print(poscar)
        
"""     
PTO
   1.00000000000000
     5.2139743000000003    0.0000000000000000    0.0000000000000000
     0.0000000000000000    5.2139743000000003    0.0000000000000000
     0.0001710487000000    0.0116047451000000    9.2051304513638925
   Pb   Ti   O
     4     4    12
Direct
  0.0000183322876524  0.0382743453515251  0.5420420037768250
  0.0000188922876490  0.0382732543515236  0.0420413847768272
  0.5000180002876463  0.5382749303515203  0.5420428407768267
  0.5000167752876480  0.5382781263515218  0.0420451497768296
  0.5000058030327477  0.0377687830048129  0.2817945210279044
  0.5000066810327493  0.0377701770048146  0.7817948410279095
  0.0000060240327481  0.5377662730048114  0.2817947600279064
  0.0000074810327514  0.5377689600048099  0.7817955580279085
  0.7500871089819571  0.2873138064180057  0.2252820891304051
  0.7500845659819563  0.2873158234180056  0.7252818861304061
  0.2500835579819581  0.7873174654180065  0.2252821191304076
  0.2500868829819596  0.7873131824180039  0.7252818261304083
  0.2499013360190072  0.2873126384180014  0.2252814321304112
  0.2498971690190070  0.2873079844180000  0.7252826031304050
  0.7498980880190084  0.7873094094180008  0.2252828361304065
  0.7499005690190111  0.7873110184180024  0.7252811691304072
  0.4999894186786342  0.0379941448076628  0.4755976959344522
  0.4999899506786356  0.0379917218076577  0.9755985839344474
  0.9999922336786361  0.5379915178076622  0.4755983529344462
  0.9999910766786363  0.5379950298076615  0.9755981059344505

"""