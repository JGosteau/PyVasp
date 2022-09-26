#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 08:01:37 2019

@author: gosteau
"""
from lib_display_BS import *
#from lib_class_Rashba import *


path = "/home/gosteau/Documents/These/2eme année/PTO/resultat/Phases/Amm2/outputs/BS_calc/2.00/bis/SOC/ST/"


class ST:
        
    def __init__(self, path, band, atom = 'tot', orb = 'tot', K=['kx','ky'], reciprocal = True, **kwargs):
        """
         - path --> chemin vers le dossier "BANDS"
         - band --> bande à traiter
         - atom = 'tot' --> projection sur l'atome atom
         - orb = 'tot'  --> projection sur les orbitales
        """
        path = clear_path(path)
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
    
    def get_k(self,kname):
        if kname == 'kx' :
            return self.kx
        elif kname == 'ky' :
            return self.ky
        elif kname == 'kz' :
            return self.kz
        
    def get_s(self,sname):
        if sname == 'sx' :
            return self.SX
        elif sname == 'sy' :
            return self.SY
        elif sname == 'sz' :
            return self.SZ
    
    def plot_ST(self,Ef = 0,K = [0,1], S = [0,1], set_color = True,
                transform_mat_k= np.eye(3),transform_mat_s= np.eye(3),
                N_interp = 0, ax = None, space = 1, **kwargs):
        if ax == None :
            ax = plt.gca()
            
            
            
        K_array = np.array([self.kx,self.ky,self.kz])
        nK = np.dot(transform_mat_k, K_array)
        KX = nK[K[0]]
        KY = nK[K[1]]
        S_array = np.array([self.SX,self.SY,self.SZ])
        nS = np.dot(transform_mat_s, S_array)
        #print(transform_mat_s)
        SX = nS[S[0]]/self.norm
        #print(self.SX - SX)
        SY = nS[S[1]]/self.norm
        #print(self.SY - SY)
        #print(S,[SX,SY])
        """
        KX = np.dot(transform_mat,self.get_k(K[0]))
        KY = np.dot(transform_mat,self.get_k(K[1]))
        
        SX = np.dot(transform_mat,self.get_s(S[0])/self.norm)
        SY = np.dot(transform_mat,self.get_s(S[1])/self.norm)
        """
        E = self.E-Ef
        if N_interp == 0 :
            if set_color == True :
                l = ax.quiver(KX[::space],KY[::space],SX[::space],SY[::space], E[::space], 
                      scale_units='xy',angles = 'xy',**kwargs)
            else :
                l = ax.quiver(KX[::space],KY[::space],SX[::space],SY[::space],  
                      scale_units='xy',angles = 'xy',**kwargs)
        
        else : 
            X = np.linspace(np.min(KX), np.max(KX), N_interp)
            Y = np.linspace(np.min(KY), np.max(KY), N_interp)
            x,y  = np.meshgrid(X, Y)
            
            sx = griddata((KX,KY), SX, (x,y))
            sz = griddata((KX,KY), SY, (x,y))
            E = griddata((KX,KY), E, (x,y))
            if set_color == True :
                l = ax.quiver(x[::space,::space],y[::space,::space], 
                          sx[::space,::space], sz[::space,::space], 
                          E[::space,::space], 
                          scale_units='xy',angles = 'xy',**kwargs)
            else :
                l = ax.quiver(x[::space,::space],y[::space,::space], 
                          sx[::space,::space], sz[::space,::space], 
                          scale_units='xy',angles = 'xy',**kwargs)
        return l
    
    def plot_ST_Ecut(self, Ecut, Ef = 0, K = [0,1], S = [0,1] ,ax= None, N_interp = 50, quiver_space = 1, args_cs = {}, 
                   decimals = 4,
                   transform_mat_k = np.eye(3), transform_mat_s = np.eye(3),
                   num_bz = [1,1], interp = False,interp_contour = False, **kwargs):
        Ecut = float(Ecut)
        if 'color' in kwargs.keys() :
            args_cs['colors'] = kwargs.get('color', None)
        if ax == None :
            ax = plt.gca()
        
        K_array = np.array([self.kx,self.ky,self.kz])
        nK = np.dot(transform_mat_k, K_array)
        
        KX = nK[K[0]]
        KY = nK[K[1]]
        
        S_array = np.array([self.SX,self.SY,self.SZ])
        nS = np.dot(transform_mat_s, K_array)
        
        SX = S_array[S[0]]/self.norm
        SY = S_array[S[1]]/self.norm
        
        #KX = np.dot(transform_mat,self.get_k(K[0]))
        #KY = np.dot(transform_mat,self.get_k(K[1]))
        
        #SX = np.dot(transform_mat,self.get_s(S[0])/self.norm)
        #SY = np.dot(transform_mat,self.get_s(S[1])/self.norm)
        
        X = np.linspace(np.min(KX), np.max(KX), N_interp)
        Y = np.linspace(np.min(KY), np.max(KY), N_interp)
        x,y  = np.meshgrid(X, Y)
        
        def get_pos_quiver(points):
            index = []
            kx = np.round(y,decimals)
            ky = np.round(x,decimals)
            #print(kx)
            #print(ky)
            for pt in np.round(points, decimals):
                #print(pt[0])
                i_kx = np.where(pt[0] == kx)[0][0]
                i_ky = np.where(pt[1] == ky)[0][0]
                #i_kx = np.where(kx == pt[0])
                #i_ky = np.where(ky == pt[1])
                #print(pt, i_kx,i_ky)
                for i in i_kx[0] :
                    if i in i_ky[0]:
                        index.append(i)
            return np.int64(index)
        
        
        E = griddata((KX,KY), self.E, (x,y))
        
        CS = ax.contour(x,y,E-Ef,Ecut, **args_cs)
        
        paths = CS.collections[0].get_paths()
        verts = [xx.vertices for xx in paths]
        points = np.concatenate(verts)
        """
        new_index = get_pos_quiver(points)
        
        kx = points[:,0]
        ky = points[:,1]
        new_sx = self.SX[new_index]
        new_sy = self.SY[new_index]
        #print(kx,ky)"""
        
        kx = points[:,0]
        ky = points[:,1]
        sx = griddata((KX,KY), SX, (kx,ky))
        sy = griddata((KX,KY), SY, (kx,ky))
        
        ax.quiver(kx,ky,sx,sy, **kwargs)
   
    def plot_ST_sz(self,Ef = 0,K = [0,1], S = [0,1,2],N_interp = 0, ax = None, space = 1, 
                   transform_mat_k= np.eye(3),transform_mat_s= np.eye(3),
                   **kwargs):
        if ax == None :
            ax = plt.gca()
        
        K_array = np.array([self.kx,self.ky,self.kz])
        nK = np.dot(transform_mat_k, K_array)
        
        KX = nK[K[0]]
        KY = nK[K[1]]
        
        S_array = np.array([self.SX,self.SY,self.SZ])
        nS = np.dot(transform_mat_s, S_array)
        
        SX = nS[S[0]]/self.norm
        SY = nS[S[1]]/self.norm
        SZ = nS[S[2]]/self.norm
        
        """
        KX = np.dot(transform_mat,self.get_k(K[0]))
        KY = np.dot(transform_mat,self.get_k(K[1]))
        
        SX = np.dot(transform_mat,self.get_s(S[0])/self.norm)
        SY = np.dot(transform_mat,self.get_s(S[1])/self.norm)
        SZ = np.dot(transform_mat,self.get_s(S[2])/self.norm)
        """
        E = self.E-Ef
        if N_interp == 0 :
            l = ax.quiver(KX[::space],KY[::space],SX[::space],SY[::space], SZ[::space], 
                      scale_units='xy',angles = 'xy',**kwargs)
        
        else : 
            X = np.linspace(np.min(KX), np.max(KX), N_interp)
            Y = np.linspace(np.min(KY), np.max(KY), N_interp)
            x,y  = np.meshgrid(X, Y)
            
            sx = griddata((KX,KY), SX, (x,y))
            sz = griddata((KX,KY), SY, (x,y))
            #E = griddata((KX,KY), E, (x,y))
            sz = griddata((KX,KY), SZ, (x,y))
        
            l = ax.quiver(x[::space,::space],y[::space,::space], 
                      sx[::space,::space], sz[::space,::space], 
                      sz[::space,::space], 
                      scale_units='xy',angles = 'xy',**kwargs)
        return l
    
   
    def plot_ST_E(self,Ef = 0,K = [0,1], S = [0,1,2], orbs = ['tot'], atoms = [0],N_interp = 0, ax = None, space = 1, 
                   transform_mat_k= np.eye(3),transform_mat_s= np.eye(3),
                   **kwargs):
        if ax == None :
            ax = plt.gca()
        
        K_array = np.array([self.kx,self.ky,self.kz])
        nK = np.dot(transform_mat_k, K_array)
        
        KX = nK[K[0]]
        KY = nK[K[1]]
        
        S_array = np.array([self.SX,self.SY,self.SZ])
        nS = np.dot(transform_mat_s, S_array)
        
        SX = nS[S[0]]/self.norm
        SY = nS[S[1]]/self.norm
        SZ = nS[S[2]]/self.norm
        
        orb = []
        for m in orbs :
            wh = np.where(self.contrib.orbs == m)[0]
            if len(wh) == 1 :
                orb.append(wh[0])
        
        tmp = np.zeros((len(self.contrib.orbs), len(KX)))
        for at in atoms :
            tmp += self.contrib.contrib[at]
        
        contrib = np.linalg.norm(tmp[orb], axis = 0)
        
        
        E = self.E-Ef
        if N_interp == 0 :
            C = contrib
            """
            l = ax.quiver(KX[::space],KY[::space],SX[::space],SY[::space], C[::space], 
                      scale_units='xy',angles = 'xy',**kwargs)"""
            l = ax.pcolor(KX[::space],KY[::space], E[::space] ,**kwargs)
            
        
        else : 
            X = np.linspace(np.min(KX), np.max(KX), N_interp)
            Y = np.linspace(np.min(KY), np.max(KY), N_interp)
            x,y  = np.meshgrid(X, Y)
            
            sx = griddata((KX,KY), SX, (x,y))
            sy = griddata((KX,KY), SY, (x,y))
            sz = griddata((KX,KY), SZ, (x,y))
            e = griddata((KX,KY), E, (x,y))
            """
            l = ax.quiver(x[::space,::space],y[::space,::space], 
                      sx[::space,::space], sy[::space,::space], 
                      C[::space,::space], 
                      scale_units='xy',angles = 'xy',**kwargs)"""
            #l = ax.pcolor(x[::space],y[::space], e[::space] ,**kwargs)
            l = ax.imshow(e, interpolation = 'nearest', extent = [np.min(x), np.max(x), np.min(y), np.max(y)],**kwargs)
        return l 
    
    def plot_ST_contrib(self,Ef = 0,K = [0,1], S = [0,1,2], orbs = ['tot'], atoms = [0],N_interp = 0, ax = None, space = 1, 
                   transform_mat_k= np.eye(3),transform_mat_s= np.eye(3),
                   **kwargs):
        if ax == None :
            ax = plt.gca()
        
        K_array = np.array([self.kx,self.ky,self.kz])
        nK = np.dot(transform_mat_k, K_array)
        
        KX = nK[K[0]]
        KY = nK[K[1]]
        
        S_array = np.array([self.SX,self.SY,self.SZ])
        nS = np.dot(transform_mat_s, S_array)
        
        SX = nS[S[0]]/self.norm
        SY = nS[S[1]]/self.norm
        SZ = nS[S[2]]/self.norm
        
        orb = []
        for m in orbs :
            wh = np.where(self.contrib.orbs == m)[0]
            if len(wh) == 1 :
                orb.append(wh[0])
        
        tmp = np.zeros((len(self.contrib.orbs), len(KX)))
        for at in atoms :
            tmp += self.contrib.contrib[at]
        
        contrib = np.linalg.norm(tmp[orb], axis = 0)
        
        
        E = self.E-Ef
        if N_interp == 0 :
            C = contrib
            """
            l = ax.quiver(KX[::space],KY[::space],SX[::space],SY[::space], C[::space], 
                      scale_units='xy',angles = 'xy',**kwargs)"""
            l = ax.pcolor(KX[::space],KY[::space], C[::space] ,**kwargs)
            
        
        else : 
            X = np.linspace(np.min(KX), np.max(KX), N_interp)
            Y = np.linspace(np.min(KY), np.max(KY), N_interp)
            x,y  = np.meshgrid(X, Y)
            
            C = griddata((KX,KY), contrib, (x,y))
            sx = griddata((KX,KY), SX, (x,y))
            sy = griddata((KX,KY), SY, (x,y))
            sz = griddata((KX,KY), SZ, (x,y))
            e = griddata((KX,KY), E, (x,y))
            """
            l = ax.quiver(x[::space,::space],y[::space,::space], 
                      sx[::space,::space], sy[::space,::space], 
                      C[::space,::space], 
                      scale_units='xy',angles = 'xy',**kwargs)"""
            l = ax.pcolor(x[::space],y[::space], C[::space] ,**kwargs)
        return l
  
