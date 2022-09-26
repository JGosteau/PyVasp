#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 09:52:02 2019

@author: gosteau
"""
from lib_class_bs import *
path = "/home/gosteau/Documents/These/1ere année/STO/STO_tetra/3.905/SOC/"
class SpT():
    def __init__(self, path, band, Ef = None, E = 0, orbs_contrib = None, ST = False, interp = None, npts = 500):
       
        if Ef == None:
            Ef = get_Ef(path, same_dir = True)
        self.Ef = Ef
        self.E = E
        self.band_sx = CONTRIB(path, band, spin = "_sx")
        self.band_sy = CONTRIB(path, band, spin = "_sy")
        self.interp = interp
        self.npts = npts
        self.ST = ST
        self.kpoints = KPOINTS(path, from_procar = True)
        
    
    def get_E_cut(self):
        if self.interp != None :
            kxmin = self.kpoints['kx'].min()
            kxmax = self.kpoints['kx'].max()
            kymin = self.kpoints['ky'].min()
            kymax = self.kpoints['ky'].max()
            kx = np.linspace(kxmin,kxmax,self.npts)
            ky = np.linspace(kymin,kymax,self.npts)
            KX,KY = np.meshgrid(kx,ky)
            z = griddata((self.kpoints['kx'], self.kpoints['ky']), self.band_sx.energy-self.Ef, (KX,KY), method = self.interp)
        else :
            tmp = pd.DataFrame()
            tmp['kx'] = self.kpoints['kx']
            tmp['ky'] = self.kpoints['ky']
            tmp['E'] = self.band_sx.energy
            return tmp
            hdf = tmp[['kx','ky','E']].pivot('kx','ky','E')
            X=hdf.columns.values
            Y=hdf.index.values
            z=hdf.values-Ef
            KX,KY = np.meshgrid(X, Y)
        fig = plt.figure('contour_tmp', figsize = (1,1))
        CS = plt.contour(KX, KY, z, [self.E], linewidths=0.5)
        paths = CS.collections[0].get_paths()
        verts = [xx.vertices for xx in paths]
        points = np.concatenate(verts)
        
