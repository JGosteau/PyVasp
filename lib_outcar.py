#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 09:18:23 2019

@author: gosteau
"""


import numpy as np

class OUTCAR():
    
    
    def __init__(self, path):
        cpath = path
        self.path = cpath
        self.read_inputs()
        
        
    def read_inputs(self):
        self.inputs = {}
        def count(string, c):
            nc = 0
            for i in string :
                if i == c:
                    nc += 1
            return nc
        outcar = open(self.path+'OUTCAR','r')
        
        for line in outcar :
            if '=' in line :
                Ar = line.split('=')
                if len(Ar) > 0:
                    for i in range(len(Ar)-1) :
                        
                        try :
                            inp = Ar[i].split()[-1]
                            out = Ar[i+1].split()[0]
                            try :
                                try :
                                    out = int(out)
                                except :
                                    out = float(out)
                            except :
                                None
                            if out == 'F' :
                                out = False
                            elif out == 'T' :
                                out = True
                            self.inputs[inp] = out
                        except :
                            None
                            
        def get_pola(self,line):
            if "p[ion]" in line :
                pion = np.array(line.split('(')[1].split(')')[0].split())
                pion = 
        
        def read_outputs(self):
            outcar = open(self.path+'OUTCAR','r')
            self.outputs = {}
            if self.inputs['LCALCPOL'] == True :
                self.POLA = {}
            for line in outcar :
                    if "p[ion]" in line :
                        
        
        outcar.close()
        


path = "/home/gosteau/Documents/These/2eme année/PTO/Structures/ortho/PBESol/1.80/BS/"