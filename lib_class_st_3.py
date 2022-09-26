#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 08:01:37 2019

@author: gosteau
"""
from lib_class_bs3 import *
#from lib_class_Rashba import *


class ST:
        
    def __init__(self, path, reciprocal = True, **kwargs):
        """
         - path --> chemin vers le dossier "BANDS"
         - band --> bande à traiter
         - atom = 'tot' --> projection sur l'atome atom
         - orb = 'tot'  --> projection sur les orbitales
        """
        self.bs = BANDS(path, **kwargs)
        
        
        """
        self.path = path
        kpoints = KPOINTS(path, **kwargs)
        kpoints.data = kpoints.data.sort_values(['kx','ky','kz'])
        
        
        
        
        self.ix = kpoints.data.index
        self.E = np.array(get_band_info(path)[0][band-1])[self.ix]
        self.contrib = CONTRIB(path, band)
        self.SX = np.array(CONTRIB(path, band, spin = '_sx')[atom][orb])[self.ix]
        self.SY = np.array(CONTRIB(path, band, spin = '_sy')[atom][orb])[self.ix]
        self.SZ = np.array(CONTRIB(path, band, spin = '_sz')[atom][orb])[self.ix]
        self.norm = np.sqrt(self.SX**2+self.SY**2+self.SZ**2)
        #ix = kpoints.data.sort_values(['kz','kx']).index
        
        if reciprocal == True :            
            self.kx = []
            self.ky = []
            self.kz = []
            for i in kpoints.data.index :
                k = kpoints['kx'][i] * kpoints.u + kpoints['ky'][i] * kpoints.v + kpoints['kz'][i] * kpoints.w
                self.kx.append(k[0])
                self.ky.append(k[1])
                self.kz.append(k[2])
            self.kx = np.array(self.kx)
            self.ky = np.array(self.ky)
            self.kz = np.array(self.kz)
        else :
            self.kx = np.array(kpoints['kx'])
            self.ky = np.array(kpoints['ky'])
            self.kz = np.array(kpoints['kz'])
        """