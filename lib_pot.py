#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 17 10:02:04 2019

@author: gosteau
"""

from lib_gen import *
from lib_poscar import *


class LOCPOT:
    def get_raw(self):
        n_pt = np.int0(grep('NGXF',self.path+'OUTCAR')[0][0].split()[-5::2])
        nx,ny,nz = n_pt
        pot=[]
        PF = grep([' '+str(nx),' '+str(ny),' '+str(nz)],self.path + self.filename)[1]
        f = open(self.path+self.filename,'r')
        count = 1
        try :
            PF[1] = PF[0]+np.ceil(nx*ny*nz/5.)
        except :
            PF.append(PF[0]+np.ceil(nx*ny*nz/5.))
        for l in f :
            if count >= PF[1]+1 :
                break
            elif count > PF[0] :
                for w in l.split() :
                    pot.append(float(w))
            count += 1
        locpot = []
        f.close()
        pot = np.array(pot)
        return pot, n_pt
    
    def xy_average(self, verb = False):
        xy_pot = np.zeros(self.nz)
        moy = self.nx*self.ny
        for z in range(self.nz):
            for i in range(moy*z,moy*(z+1)):
                xy_pot[z] += self.data[i]/moy
        return xy_pot
    
    def save(self):
        np.savez_compressed(self.path+self.filename.lower(),**self.__dict__)
        
    def load(self):
        loaded = np.load(self.path+self.filename.lower()+'.npz')
        for key in loaded.__dict__['files']:
            self.__dict__[key] = loaded[key]
            try :
                len(self.__dict__[key])
            except:
                if 'int' in self.__dict__[key].dtype.name :
                    self.__dict__[key] = int(self.__dict__[key])
                elif 'float' in self.__dict__[key].dtype.name :
                    self.__dict__[key] = float(self.__dict__[key])
                elif 'string' in self.__dict__[key].dtype.name :
                    self.__dict__[key] = str(self.__dict__[key])
        
    
    def __init__(self, path, filename = 'LOCPOT', LOAD = True, poscar = None):
        self.path = clear_path(path,space=False)
        self.filename = filename
        self.data = None 
        self.nx = None
        self.ny = None
        self.nz = None
        self.X = None
        self.Y = None
        self.Z = None
        self.diag_lattice = None
        
        if LOAD == True and os.path.isfile(self.path+self.filename.lower()+'.npz') == True:
            print ('LOADING')
            self.load()
        else :
            self.data, npt = self.get_raw()
            if poscar == None:
                poscar = POSCAR(self.path)
            self.diag_lattice = poscar.diag_lattice
            self.nx, self.ny, self.nz = npt
            self.X = np.array(range(self.nx)*self.ny*self.nz)
            
            Y = [[i]*self.nx for i in range(self.ny)]
            self.Y = []
            for i in Y:
                self.Y += i
            self.Y = np.array(self.Y*self.nz)
            
            Z = [[i]*self.nx*self.ny for i in range(self.nz)]
            self.Z = []
            for i in Z:
                self.Z += i
            self.Z = np.array(self.Z)
            self.save()
            
    
            
    def __getitem__(self,item):
        if item == 0 or item == 'x' or item == 'X':
            return self.X
        elif item == 1 or item == 'y' or item == 'Y':
            return self.Y
        elif item == 2 or item == 'z' or item == 'Z':
            return self.Z
        else:
            return self.data

    def search_closest_index(self,z,iprec,z_moy,dz_moy):
        if z < z_moy - dz_moy/2 :
            return iprec - 1
        elif z > z_moy + dz_moy/2 :
            return iprec + 1
        else :
            return iprec

    def moy_z_pot(self,xy_pot, z_moy, dz):
        moy = 0.
        count = 0
        
        try :
            if z_moy+dz+1 > self.nz :
                R = list(range(z_moy-dz, self.nz))+list(range(0,(z_moy + dz + 1)%self.nz))
            else :
                R = list(range(z_moy-dz, z_moy + dz + 1))
            
            for z in R :
                moy += xy_pot[z]
                count += 1
            return moy/count
        except :
            print( z_moy-dz, z_moy+dz+1,(z_moy-dz)%self.nz, (z_moy + dz+1)%self.nz, self.nz)
            
            if z_moy+dz+1 > self.nz :
                R = list(range(z_moy-dz, self.nz))+list(range(0,(z_moy + dz + 1)%self.nz))
            else :
                R = list(range((z_moy-dz), (z_moy + dz+1)))
            print( R)
            print( xy_pot[R])
            raise TypeError
            
        

    def macro_average(self, poscar = None, verb = False, scale = 1, XY_pot = None, dist_planes = None) :
        if poscar == None:
            poscar = POSCAR(self.path)
        ax,ay,az = poscar.diag_lattice
        if type(dist_planes) == type(None):
            dist_planes = poscar.get_dist_plane()
        if type(XY_pot) == type(None):
            XY_pot = self.xy_average()
        Z_moy = np.int0(dist_planes['z'] * self.nz)
        DZ_moy = np.int0(dist_planes['dz'] * self.nz)
        Z_pot_moy = np.zeros(self.nz)
        index = 0
        for z in range(0,self.nz): 
            if z < Z_moy[0] - scale * DZ_moy[0] or z > Z_moy[-1] + scale * DZ_moy[-1] :
                Z_pot_moy[z] = XY_pot[z]
            else :
                index = self.search_closest_index(z, index, Z_moy[index], scale * DZ_moy[index])%len(DZ_moy)
                Z_pot_moy[z] = self.moy_z_pot(XY_pot, z, scale * DZ_moy[index])
                    
        return Z_pot_moy
    
    def get_Efield(self, POT, zi, zf, verb = True, display = False, color = 'black', dz_scale = 2 ,display_text = True ,display_text_pos = 'right'):
        fit = lambda x,a,b : a*x+b
        Z = np.arange(zi,zf+1)
        nz = self.nz
        az = self.diag_lattice[-1]
        popt, _ = curve_fit(fit,Z,POT[zi:zf+1])
        if verb :
            print( "Efield = %e eV/A" %(-popt[0]*nz/az))
        if display :
            text = "%1.2e eV/A" %(-popt[0]*nz/az)
            dz = (zf-zi)*dz_scale
            Z2 = np.linspace(zi-dz,zf+dz, 1000)
            plt.plot(Z2,fit(Z2,*popt), color = color, ls = '--')
            if display_text :
                if display_text_pos == 'left':
                    plt.text(Z[0], fit(Z[0],*popt), text, color = color,bbox=dict(boxstyle="square",fc = 'white'))
                else :
                    plt.text(Z[-1], fit(Z[-1],*popt), text, color = color,bbox=dict(boxstyle="square",fc = 'white'))
        return -popt[0]*nz/az

def PLOT_LOC_vs_PLANES(path = None, locpot = None, poscar = None, z_range = None, dist_planes = None, ax = None, ylim = None,**kwargs):
    if locpot == None:
        locpot = LOCPOT(path)
    if poscar == None:
        poscar = POSCAR(path)
    if type(dist_planes) == type(None):
        dist_planes = poscar.get_dist_plane(d = 1.0)
    xy_pot = locpot.xy_average()
    macro = locpot.macro_average(poscar = poscar,XY_pot = xy_pot, dist_planes = dist_planes)
    if z_range == None :
        z_range = range(locpot.nz)
    if ax == None :
        ax = plt.subplot()
    plt.plot(xy_pot, *kwargs)
    plt.plot(macro, *kwargs)
    if ylim == None :
        plt.ylim(*plt.ylim())
    else :
        plt.ylim(*ylim)
    plt.xlim(0,locpot.nz)
    Z_planes = np.int0(dist_planes['z'] * locpot.nz)
    plt.vlines(np.int0(dist_planes['z'] * locpot.nz), *plt.ylim(),color = 'lightgrey')
    #ax2 = ax.twiny()
    #ax2.get_xaxis().set_visible(False)
    #ax.set_xticks(Z_planes[::2])
    #ax.set_xticklabels(range(0,len(Z_planes),2))
    #ax2.set_xbound(0,locpot.nz)
    return macro
    
    
    
    
        
    
    
        
