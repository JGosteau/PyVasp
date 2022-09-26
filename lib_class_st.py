#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 08:01:37 2019

@author: gosteau
"""
from lib_display_BS import *
from lib_class_Rashba import *

mpl.rcParams['text.usetex'] = False
class ST:
    def __init__(self, path, band, atom = 'tot', orb = 'tot', K=['kx','ky'], **kwargs):
        """
         - path --> chemin vers le dossier "BANDS"
         - band --> bande à traiter
         - atom = 'tot' --> projection sur l'atome atom
         - orb = 'tot'  --> projection sur les orbitales
        """
        path = clear_path(path)
        self.path = path
        E = get_band_info(path)[0][band-1]
        kpoints = KPOINTS(path, **kwargs)
        SX = CONTRIB(path, band, spin = '_sx')[atom][orb]
        SY = CONTRIB(path, band, spin = '_sy')[atom][orb]
        SZ = CONTRIB(path, band, spin = '_sz')[atom][orb]
        self.ix = kpoints.data.sort_values([K[1],K[0]]).index
        #ix = kpoints.data.sort_values(['kz','kx']).index
        self.kx, self.ky, self.kz = np.meshgrid(np.unique(kpoints.loc[self.ix][K[0]]),np.unique(kpoints.loc[self.ix][K[1]],np.unique(kpoints.loc[self.ix][K[2]])))
        #self.kx, self.ky = np.meshgrid(np.unique(kpoints.loc[ix]['kx']),np.unique(kpoints.loc[ix]['kz']))
        self.lx = len(np.unique(kpoints.loc[self.ix][K[0]]))
        self.ly = len(np.unique(kpoints.loc[self.ix][K[1]]))
        #lx = len(np.unique(kpoints.loc[ix]['kx']))
        #ly = len(np.unique(kpoints.loc[ix]['kz']))
        self.sx = SX[self.ix].reshape(self.ly,self.lx)
        self.sy = SY[self.ix].reshape(self.ly,self.lx)
        self.sz = SZ[self.ix].reshape(self.ly,self.lx)
        self.norm = np.sqrt(self.sx**2+self.sy**2+self.sz**2)
        self.E = E[self.ix].reshape(self.ly,self.lx)
        #self.get_reciprocal_lattice_vectors()
    
    
    def get_kx_and_ky(self):
        kpoints = self.KPOINTS(path, fileformat = 'PROCAR')
        self.ix = kpoints.data.sort_values()
    
    def get_reciprocal_lattice_vectors(self):
        file = path + 'OUTCAR'
        n = grep('reciprocal lattice vectors', file)[1][0]
        n = get_row(n+1, file).split()
        self.b1 = np.float64(get_row(n+1, file).split())[3:]
        self.b2 = np.float64(get_row(n+2, file).split())[3:]
        self.b3 = np.float64(get_row(n+3, file).split())[3:]
        
        
        
    
    def plot_ST(self,Ef = 0, ax = None, quiver_space = 1, S=['sx','sy'],**kwargs):
        if ax == None :
            ax = plt
        if S[0] == 'sx' :
            S1 = self.sx
        if S[0] == 'sy' :
            S1 = self.sy
        if S[0] == 'sz' :
            S1 = self.sz
        if S[1] == 'sx' :
            S2 = self.sx
        if S[1] == 'sy' :
            S2 = self.sy
        if S[1] == 'sz' :
            S2 = self.sz
        l = ax.quiver(self.kx[::quiver_space,::quiver_space],self.ky[::quiver_space,::quiver_space],
           (S1/self.norm)[::quiver_space,::quiver_space],(S2/self.norm)[::quiver_space,::quiver_space], 
           self.E[::quiver_space,::quiver_space]-Ef,
           scale_units='xy',angles = 'xy', 
           **kwargs)
        return l
    
        
    
    
    def plot_ST_Ecut(self, E, Ef = 0, ax = None, quiver_space = 1, args_cs = {},S=['sx','sy'], num_bz = [1,1], **kwargs):
        if 'color' in kwargs.keys() :
            args_cs['colors'] = kwargs.get('color', None)
        #colors = args_cs.get('colors', None)
        #color = kwargs.get('color', colors)
        #kwargs['color'] = color
        
        if S[0] == 'sx' :
            S1 = self.sx
        if S[0] == 'sy' :
            S1 = self.sy
        if S[0] == 'sz' :
            S1 = self.sz
        if S[1] == 'sx' :
            S2 = self.sx
        if S[1] == 'sy' :
            S2 = self.sy
        if S[1] == 'sz' :
            S2 = self.sz
        def get_pos_quiver(points):
            index = []
            decimals = 2
            kx = np.round(self.kx,decimals)
            ky = np.round(self.ky,decimals)
            #print(kx)
            #print(ky)
            for pt in np.round(points, decimals):
                #print(pt[0])
                #i_kx = np.where(pt[0] == kx)[0][0]
                #i_ky = np.where(pt[1] == ky)[0][0]
                i_kx = np.searchsorted(kx[0], pt[0])
                i_ky = np.searchsorted(ky[:,0], pt[1])
                
                index.append([i_kx,i_ky])
            return np.int64(index)
    
        if ax == None :
            ax = plt
        CS = []
        
        for i in np.arange(num_bz[0]):
            for j in np.arange(num_bz[1]):
                CS.append(ax.contour(self.kx+i,self.ky+j,self.E-Ef, E, **args_cs))
        CS = CS[0]
        
        paths = CS.collections[0].get_paths()
        verts = [xx.vertices for xx in paths]
        points = np.concatenate(verts)
        ix = get_pos_quiver(points)
        
        new_kx_ix = ix[:,0]
        new_ky_ix = ix[:,1]
        new_norm = self.norm[new_ky_ix,new_kx_ix]
        new_Sx = S1[new_ky_ix,new_kx_ix]/new_norm
        new_Sy = S2[new_ky_ix,new_kx_ix]/new_norm
        for i in np.arange(num_bz[0]):
            for j in np.arange(num_bz[1]):
                
                ax.quiver(points[::quiver_space,0]+i,points[::quiver_space,1]+j,new_Sx[::quiver_space],new_Sy[::quiver_space], **kwargs)
        
        
    def plot_ST_contrib_z(self, Ef = 0, ax = None, quiver_space = 1, args_cs = {},S=['sx','sy','sz'], **kwargs):
        if ax == None :
            ax = plt
        #colors = args_cs.get('colors', None)
        #color = kwargs.get('color', colors)
        #kwargs['color'] = color
        
        if S[0] == 'sx' :
            S1 = self.sx
        if S[0] == 'sy' :
            S1 = self.sy
        if S[0] == 'sz' :
            S1 = self.sz
        if S[1] == 'sx' :
            S2 = self.sx
        if S[1] == 'sy' :
            S2 = self.sy
        if S[1] == 'sz' :
            S2 = self.sz
        if S[2] == 'sx' :
            S3 = self.sx
        if S[2] == 'sy' :
            S3 = self.sy
        if S[2] == 'sz' :
            S3 = self.sz
            
        l = ax.quiver(self.kx[::quiver_space,::quiver_space],self.ky[::quiver_space,::quiver_space],
           (S1/self.norm)[::quiver_space,::quiver_space],(S2/self.norm)[::quiver_space,::quiver_space], 
           (S3/self.norm)[::quiver_space,::quiver_space],
           scale_units='xy',angles = 'xy', 
           **kwargs)
        return l

class ST:
        
    def __init__(self, path, band, atom = 'tot', orb = 'tot', K=['kx','ky'], **kwargs):
        """
         - path --> chemin vers le dossier "BANDS"
         - band --> bande à traiter
         - atom = 'tot' --> projection sur l'atome atom
         - orb = 'tot'  --> projection sur les orbitales
        """
        path = clear_path(path)
        self.path = path
        E = get_band_info(path)[0][band-1]
        kpoints = KPOINTS(path, **kwargs)
        SX = CONTRIB(path, band, spin = '_sx')[atom][orb]
        SY = CONTRIB(path, band, spin = '_sy')[atom][orb]
        SZ = CONTRIB(path, band, spin = '_sz')[atom][orb]
        
        
        
        
        self.ix = kpoints.data.sort_values([K[1],K[0]]).index
        #ix = kpoints.data.sort_values(['kz','kx']).index
        self.kx, self.ky = np.meshgrid(np.unique(kpoints.loc[self.ix][K[0]]),np.unique(kpoints.loc[self.ix][K[1]]))
        #self.kx, self.ky = np.meshgrid(np.unique(kpoints.loc[ix]['kx']),np.unique(kpoints.loc[ix]['kz']))
        self.lx = len(np.unique(kpoints.loc[self.ix][K[0]]))
        self.ly = len(np.unique(kpoints.loc[self.ix][K[1]]))
        normal = False
        if normal :
        
        #lx = len(np.unique(kpoints.loc[ix]['kx']))
        #ly = len(np.unique(kpoints.loc[ix]['kz']))
            self.sx = SX[self.ix].reshape(self.ly,self.lx)
            self.sy = SY[self.ix].reshape(self.ly,self.lx)
            self.sz = SZ[self.ix].reshape(self.ly,self.lx)
            self.norm = np.sqrt(self.sx**2+self.sy**2+self.sz**2)
            self.E = E[self.ix].reshape(self.ly,self.lx)
            self.norm = np.sqrt(self.sx**2+self.sy**2+self.sz**2)
            #self.get_reciprocal_lattice_vectors()
        else : 
            self.kx, self.ky = np.meshgrid(np.unique(kpoints.loc[self.ix][K[0]]),np.unique(kpoints.loc[self.ix][K[1]]))
            self.sx = np.zeros((self.ly, self.lx))    
            self.sy = np.zeros((self.ly, self.lx))    
            self.sz = np.zeros((self.ly, self.lx))   
            self.E = np.zeros((self.ly, self.lx))      
        
            import copy
            kx_table = copy.deepcopy(self.kx[0])
            ky_table = copy.deepcopy(self.ky[:,0])
            for j in range(len(self.ky[0])):
                kx = kx_table[j]
                for i in range(len(self.ky)):
                    ky = ky_table[i]
                    new_kpt = kpoints.loc[kpoints[K[0]] == kx].loc[kpoints[K[1]] == ky]
                    cond = len(new_kpt)
                    #print(i,j,kx,ky,cond)
                    if cond == 0:
                        self.kx[i][j] = np.nan
                        self.ky[i][j] = np.nan 
                        self.sx[i][j] = np.nan
                        self.sy[i][j] = np.nan
                        self.sz[i][j] = np.nan
                    else :
                        ix = new_kpt.index[0]
                        self.sx[i][j] = SX[ix]
                        self.sy[i][j] = SY[ix]
                        self.sz[i][j] = SZ[ix]
                        self.E[i][j] = E[ix]
            self.norm = np.sqrt(self.sx**2+self.sy**2+self.sz**2)
                    
    
        
    
    def get_kx_and_ky(self):
        kpoints = self.KPOINTS(path, fileformat = 'PROCAR')
        self.ix = kpoints.data.sort_values()
    
    def get_reciprocal_lattice_vectors(self):
        file = path + 'OUTCAR'
        n = grep('reciprocal lattice vectors', file)[1][0]
        n = get_row(n+1, file).split()
        self.b1 = np.float64(get_row(n+1, file).split())[3:]
        self.b2 = np.float64(get_row(n+2, file).split())[3:]
        self.b3 = np.float64(get_row(n+3, file).split())[3:]
        
        
        
    
    def plot_ST(self,Ef = 0, ax = None, quiver_space = 1, S=['sx','sy'],**kwargs):
        if ax == None :
            ax = plt
        if S[0] == 'sx' :
            S1 = self.sx
        if S[0] == 'sy' :
            S1 = self.sy
        if S[0] == 'sz' :
            S1 = self.sz
        if S[1] == 'sx' :
            S2 = self.sx
        if S[1] == 'sy' :
            S2 = self.sy
        if S[1] == 'sz' :
            S2 = self.sz
        l = ax.quiver(self.kx[::quiver_space,::quiver_space],
                      self.ky[::quiver_space,::quiver_space],
                      (S1/self.norm)[::quiver_space,::quiver_space],
                      (S2/self.norm)[::quiver_space,::quiver_space],
                      self.E[::quiver_space,::quiver_space]-Ef,
           scale_units='xy',angles = 'xy', 
           **kwargs)
        return l
    
        
    
    
    def plot_ST_Ecut(self, E, Ef = 0, ax = None, quiver_space = 1, args_cs = {},S=['sx','sy'], num_bz = [1,1], interp = False,interp_contour = False, **kwargs):
        if 'color' in kwargs.keys() :
            args_cs['colors'] = kwargs.get('color', None)
        #colors = args_cs.get('colors', None)
        #color = kwargs.get('color', colors)
        #kwargs['color'] = color
        
        if S[0] == 'sx' :
            S1 = self.sx
        if S[0] == 'sy' :
            S1 = self.sy
        if S[0] == 'sz' :
            S1 = self.sz
        if S[1] == 'sx' :
            S2 = self.sx
        if S[1] == 'sy' :
            S2 = self.sy
        if S[1] == 'sz' :
            S2 = self.sz
        def get_pos_quiver(points):
            index = []
            decimals = 4
            kx = np.round(self.kx,decimals)
            ky = np.round(self.ky,decimals)
            #print(kx)
            #print(ky)
            for pt in np.round(points, decimals):
                #print(pt[0])
                #i_kx = np.where(pt[0] == kx)[0][0]
                #i_ky = np.where(pt[1] == ky)[0][0]
                i_kx = np.searchsorted(kx[0], pt[0])
                i_ky = np.searchsorted(ky[:,0], pt[1])
                
                index.append([i_kx,i_ky])
            return np.int64(index)
    
        if ax == None :
            ax = plt
        CS = []
        
        if interp_contour == True :
            x0 = np.min(self.kx)
            xf = np.max(self.kx)
            y0 = np.min(self.kx)
            yf = np.max(self.kx)
            n_new = 101
            new_kx, new_ky = np.meshgrid(np.linspace(x0,xf,n_new),np.linspace(y0,yf,n_new))
            kx = self.kx[np.isnan(self.kx) == False]
            ky = self.ky[np.isnan(self.ky) == False]
            old_E = self.E[np.isnan(self.kx) == False]
            new_E = griddata((kx,ky), old_E, (new_kx,new_ky))
            print(new_E)
            CS=ax.contour(new_kx,new_ky,new_E-Ef, E, **args_cs)
            print(CS, CS.collections[0].get_paths())
            """
            for i in np.arange(num_bz[0]):
                    for j in np.arange(num_bz[1]):
                        CS.append(ax.contour(new_kx+i,new_ky+j,new_E-Ef, E, **args_cs))"""
        else :        
            for i in np.arange(num_bz[0]):
                for j in np.arange(num_bz[1]):
                    CS.append(ax.contour(self.kx+i,self.ky+j,self.E-Ef, E, **args_cs))
        CS = CS[0]
        
        paths = CS.collections[0].get_paths()
        verts = [xx.vertices for xx in paths]
        points = np.concatenate(verts)
        KX = points[:,0][np.isnan(points[:,1]) == False]
        KY = points[:,1][np.isnan(points[:,1]) == False]
        if interp == False :
            ix = get_pos_quiver(points)
            """ proposer interpolation """
            new_kx_ix = ix[:,0]
            new_ky_ix = ix[:,1]
            new_norm = self.norm[new_ky_ix,new_kx_ix]
            new_Sx = S1[new_ky_ix,new_kx_ix]/new_norm
            new_Sy = S2[new_ky_ix,new_kx_ix]/new_norm
        else :
            points = np.meshgrid(KX,KY)
            #return KX,KY
            tmp_sx = np.zeros((self.ly,self.lx))
            tmp_sy = np.zeros((self.ly,self.lx))
            method = 'cubic'
            sx = self.sx[np.isnan(self.sx) == False]
            sy = self.sy[np.isnan(self.sy) == False]
            kx = self.kx[np.isnan(self.kx) == False]
            ky = self.ky[np.isnan(self.ky) == False]
            new_Sx = griddata((kx,ky), sx, (KX,KY))
            new_Sy = griddata((kx,ky), sy, (KX,KY))
        for i in np.arange(num_bz[0]):
            for j in np.arange(num_bz[1]):
                
                #ax.quiver(points[::quiver_space,0]+i,points[::quiver_space,1]+j,new_Sx[::quiver_space],new_Sy[::quiver_space], **kwargs)
                ax.quiver(KX[::quiver_space]+i,KY[::quiver_space]+j,new_Sx[::quiver_space],new_Sy[::quiver_space], **kwargs)
        
        
    def plot_ST_contrib_z(self, Ef = 0, ax = None, quiver_space = 1, args_cs = {},S=['sx','sy','sz'], **kwargs):
        if ax == None :
            ax = plt
        #colors = args_cs.get('colors', None)
        #color = kwargs.get('color', colors)
        #kwargs['color'] = color
        
        if S[0] == 'sx' :
            S1 = self.sx
        if S[0] == 'sy' :
            S1 = self.sy
        if S[0] == 'sz' :
            S1 = self.sz
        if S[1] == 'sx' :
            S2 = self.sx
        if S[1] == 'sy' :
            S2 = self.sy
        if S[1] == 'sz' :
            S2 = self.sz
        if S[2] == 'sx' :
            S3 = self.sx
        if S[2] == 'sy' :
            S3 = self.sy
        if S[2] == 'sz' :
            S3 = self.sz
            
        l = ax.quiver(self.kx[::quiver_space,::quiver_space],self.ky[::quiver_space,::quiver_space],
           (S1/self.norm)[::quiver_space,::quiver_space],(S2/self.norm)[::quiver_space,::quiver_space], 
           (S3/self.norm)[::quiver_space,::quiver_space],
           scale_units='xy',angles = 'xy', 
           **kwargs)
        return l
        
FIT = {'GX' : lambda x,g, g2, E0 : 2*g*x+2*g2*x**3 + E0, 'GM' : lambda x,g, gt, E0 : 2*g*x+gt*x**3 + E0}
def rashba_p4mm(path,b1,b2,R, fit, ref_kpt = None, ax = 1, ay = 1, az= 1, display = True, fit_band = True):
    
    A = {'kx' : ax, 'ky' : ay, 'kz' : az}
    band_1 = get_band_info(path)[0][b1-1][R]
    band_2 = get_band_info(path)[0][b2-1][R]
    
    
    kpoints = KPOINTS(path,fileformat='PROCAR')
    new_kpoints = kpoints.loc[R].reset_index()
    if ref_kpt == None:
        ref_kpt = new_kpoints.loc[0]
    k_norm = np.zeros(len(R))
    for k in ['kx','ky','kz']:
        k_norm += ((new_kpoints[k]-ref_kpt[k])/A[k])**2
    k_norm = np.sqrt(k_norm)
    if display :
        ix = 0
        kpt = R[ix]
        kx = np.array(new_kpoints['kx'])[ix]
        ky = np.array(new_kpoints['ky'])[ix]
        kz = np.array(new_kpoints['kz'])[ix]
        print("fit de kpt %4d : kx = %.8f ; ky = %.8f ; kz = %.8f" %(kpt, kx,ky,kz))
        ix = -1
        kpt = R[ix]
        kx = np.array(new_kpoints['kx'])[ix]
        ky = np.array(new_kpoints['ky'])[ix]
        kz = np.array(new_kpoints['kz'])[ix]
        print("à      kpt %4d : kx = %.8f ; ky = %.8f ; kz = %.8f" %(kpt, kx,ky,kz))
    
    popt, err = curve_fit(fit, k_norm, band_2-band_1)
    err = np.sqrt(np.diag(err))/popt*100
    #err = np.sqrt(np.diag(err))*100
    #for i in range(len(popt)):
    #    err[i] /= popt[i]
    
    if fit_band == True :
        
        fit_band_2 = lambda x,a : a*x**2 + 1/2*fit(x, *popt) + band_2[0]
        fit_band_1 = lambda x,a : a*x**2 - 1/2*fit(x, *popt) + band_1[0]
        popt_a1, err_a1 = curve_fit(fit_band_1, k_norm, band_1)
        err_a1 = np.sqrt(np.diag(err_a1))/popt_a1*100
        a1 = popt_a1[0]
        ea1 = err_a1[0]
        popt_a2, err_a2 = curve_fit(fit_band_1, k_norm, band_2)
        err_a2 = np.sqrt(np.diag(err_a2))/popt_a2*100
        a2 = popt_a2[0]
        ea2 = err_a2[0]
        return popt, err, a1,ea1,a2,ea2
        
        return popt, err, 0,0,0,0
    else :
        return popt, err



def fit_to_choose(g, g2, E0, direction = 'GX'):
    A = np.array([g,g2,E0])
    ix = np.where(A != None)[0]
    if len(ix) == 0 :
        ix_VAL = {'g' : 0, 'g2' : 1,'E0' : 2}
        l = FIT[direction]
    elif len(ix) == 1 :
        if E0 != None :
            ix_VAL = {'g' : 0, 'g2' : 1}
            l = lambda x, g,g2 : FIT[direction](x,g,g2,E0)
        elif g != None:
            ix_VAL = {'g2' : 0, 'E0' : 1}
            l = lambda x, g2,E0 : FIT[direction](x,g,g2,E0)
        elif g2 != None:
            ix_VAL = {'g' : 0, 'E0' : 1}
            l = lambda x, g,E0 : FIT[direction](x,g,g2,E0)
    elif len(ix) == 2 :
        if E0 != None and g != None :
            ix_VAL = {'g2' : 0}
            l = lambda x, g2 : FIT[direction](x,g,g2,E0)
        elif E0 != None and g2 != None :
            ix_VAL = {'g' : 0}
            l = lambda x, g : FIT[direction](x,g,g2,E0)
        elif g2 != None and g1 != None :
            ix_VAL = {'E0' : 0}
            l = lambda x, E0 : FIT[direction](x,g,g2,E0)
    else :
        ix_VAL = {}
        l = lambda x, a : x*a
    return l, ix_VAL



def rashba_auto_p4mm(path, b1, b2, R_GX, R_GM,poscar = None, kpoints = None, E0 = None, ref_kpt = None, display = True, g = None, g1 = None, g2 = None):
    poscar = POSCAR(path)
    ax, ay, az = poscar.diag_lattice
    R = R_GX
    #ref_kpt = {'kx' : 0, 'ky' : 0, 'kz' : 0}
    A = np.array([ax,ay,az])/(2*np.pi)
    band_1 = get_band_info(path)[0][b1-1][R]
    band_2 = get_band_info(path)[0][b2-1][R]
    kpoints = KPOINTS(path,fileformat='PROCAR')
    new_kpoints = kpoints.loc[R].reset_index()
    if ref_kpt == None:
        ref_kpt = new_kpoints.loc[0]
    k_norm = np.zeros(len(R))
    for k in ['kx','ky','kz']:
        k_norm += ((new_kpoints[k]-ref_kpt[k]))**2
    k_norm = np.sqrt(k_norm)
    # GX :
    VAR = ['g','g2','E0']
    INI_VAR = {'g' : g,'g2' : g2,'E0' : E0}
    fit, ix_VAR = fit_to_choose(g,g2,E0)
    def search_COEFF(popt, err, ix_VAR, INI_VAR):
        COEFF = {}
        ERROR = {}
        for V in VAR :
            if V in ix_VAR.keys():
                COEFF[V] = popt[ix_VAR[V]]
                ERROR[V] = err[ix_VAR[V]]
            else :
                COEFF[V] = INI_VAR[V]
                ERROR[V] = 0
        return COEFF, ERROR
    
    popt_GX_norm, err_GX_norm, a1_GX_norm, ea1_GX_norm, a2_GX_norm, ea2_GX_norm = rashba_p4mm(path, b1,b2,R_GX, fit, ref_kpt, *A, display = display)#ax = A[0], ay = A[1], az = A[2])
    popt_GX_ua, err_GX_ua, a1_GX_ua, ea1_GX_ua, a2_GX_ua, ea2_GX_ua = rashba_p4mm(path, b1,b2,R_GX, fit, ref_kpt, display = False)
    GX_norm, GX_error_norm = search_COEFF(popt_GX_norm, err_GX_norm, ix_VAR, INI_VAR)
    GX_ua, GX_error_ua = search_COEFF(popt_GX_ua, err_GX_ua, ix_VAR, INI_VAR)
    if display :
        print("              --------------------- GX ---------------------")
        print("fit GX : dE = 2*g*kx + 2*g2*kx**3 + E0")
    
    if display :
        g = GX_norm['g']
        eg = GX_error_norm['g']
        g2 = GX_norm['g2']
        g2_norm = g2
        eg2 = GX_error_norm['g2']
        E0 = GX_norm['E0']
        eE0 = GX_error_norm['E0']
        
        
        a1 =  a1_GX_norm
        ea1 = ea1_GX_norm
        a2 =  a2_GX_norm
        ea2 = ea2_GX_norm
        print("Normalisé : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g2 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2, g,eg,g2,eg2,E0,eE0))
            
    if display :
        g = GX_ua['g']
        eg = GX_error_ua['g']
        g2 = GX_ua['g2']
        eg2 = GX_error_ua['g2']
        E0 = GX_ua['E0']
        eE0 = GX_error_ua['E0']
        g2_ua = g2
        
        a1 =  a1_GX_ua
        ea1 = ea1_GX_ua
        a2 =  a2_GX_ua
        ea2 = ea2_GX_ua
        print("ua        : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g2 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2,g,eg,g2,eg2,E0,eE0))
        print("              ----------------------------------------------")
    #plt.plot(k_norm, band_2-band_1, color = 'black')
    #plt.plot(k_norm, FIT['GX'](k_norm, *popt_GX_ua,E0),color = 'red', linestyle = ':')
    # GM :
    VAR = ['g','g2','E0']
    g = INI_VAR['g']
    E0 = INI_VAR['E0']
    
    if g1 != None :
        try :
            INI_VAR = {'g' : g,'g2' : g2_norm,'E0' : E0}
        except :
            INI_VAR = {'g' : g,'g2' : 0,'E0' : E0}
            
    else :
        INI_VAR = {'g' : g,'g2' : g1,'E0' : E0}
    fit, ix_VAR = fit_to_choose(g,g1,E0, direction = 'GM')
    popt_GM_norm, err_GM_norm, a1_GM_norm, ea1_GM_norm, a2_GM_norm, ea2_GM_norm = rashba_p4mm(path, b1,b2,R_GM, fit, ref_kpt, *A, display = display)
    GM_norm, GM_error_norm = search_COEFF(popt_GM_norm, err_GM_norm, ix_VAR, INI_VAR)
    
    if g1 != None :
        try :
            INI_VAR = {'g' : g,'g2' : g2_ua,'E0' : E0}
        except :
            INI_VAR = {'g' : g,'g2' : 0,'E0' : E0}
    fit, ix_VAR = fit_to_choose(g,g1,E0, direction = 'GM')
    popt_GM_ua, err_GM_ua, a1_GM_ua, ea1_GM_ua, a2_GM_ua, ea2_GM_ua = rashba_p4mm(path, b1,b2,R_GM, fit, ref_kpt, display = False)
    GM_ua, GM_error_ua = search_COEFF(popt_GM_ua, err_GM_ua, ix_VAR, INI_VAR)
    if display :
        print("              --------------------- GM ---------------------")
        print("fit GM : dE = 2*g*kx + 2*(g2+g3)*kx**3 + E0")
    
    GM_norm['g2'] -= GX_norm['g2']
    GM_error_norm['g2'] += GX_error_norm['g2']
    if display :
        g = GM_norm['g']
        eg = GM_error_norm['g']
        g1 = GM_norm['g2']
        if INI_VAR['g2'] == None :
            eg1 = GM_error_norm['g2']
        else :
            eg1 = 0
        E0 = GM_norm['E0']
        eE0 = GM_error_norm['E0']
        
        a1 =  a1_GM_norm
        ea1 = ea1_GM_norm
        a2 =  a2_GM_norm
        ea2 = ea2_GM_norm
        print("Normalisé : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g1 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2,g,eg,g1,eg1,E0,eE0))
    
    GM_ua['g2'] -= GX_ua['g2']
    GM_error_ua['g2'] += GX_error_ua['g2']
    if display :
        g = GM_ua['g']
        eg = GM_error_ua['g']
        g1 = GM_ua['g2']
        if INI_VAR['g2'] == None :
            eg1 = GM_error_ua['g2']
        else :
            eg1 = 0
        E0 = GM_ua['E0']
        eE0 = GM_error_ua['E0']
        
        a1 =  a1_GM_ua
        ea1 = ea1_GM_ua
        a2 =  a2_GM_ua
        ea2 = ea2_GM_ua
        print("ua        : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g1 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2,g,eg,g1,eg1,E0,eE0))
        print("              ----------------------------------------------")
    RES = {}
    RES['norm'] ={'GX' :  {'a1' : a1_GX_norm , 'a2' : a2_GX_norm, 'g': GX_norm['g'],     'g2' : GX_norm['g2']}, 
                  'GM' :  {'a1' : a1_GM_norm , 'a2' : a2_GM_norm, 'g': GM_norm['g'],     'g1' : GM_norm['g2']}}
    RES['ua'] ={'GX' :    {'a1' : a1_GX_ua   , 'a2' : a2_GX_ua,   'g': GX_ua['g'],       'g2' : GX_ua['g2']}, 
                'GM' :    {'a1' : a1_GM_ua   , 'a2' : a2_GM_ua,   'g': GM_ua['g'],       'g1' : GM_ua['g2']}}
    RES['error'] ={'GX' : {'a1' : ea1_GX_ua  , 'a2' : ea2_GX_ua,  'g': GX_error_ua['g'], 'g2' : GX_error_ua['g2']}, 
                   'GM':  {'a1' : ea1_GM_ua  , 'a2' : ea2_GM_ua,  'g': GM_error_ua['g'], 'g1' : GM_error_ua['g2']}}
    return RES
"""
def rashba_auto_p4mm(path, b1, b2, R_GX, R_GM,poscar = None, kpoints = None, E0 = None, ref_kpt = None, display = True, g = None, g1 = None, g2 = None):
    poscar = POSCAR(path)
    ax, ay, az = poscar.diag_lattice
    R = R_GX
    #ref_kpt = {'kx' : 0, 'ky' : 0, 'kz' : 0}
    A = np.array([ax,ay,az])/(2*np.pi)
    band_1 = get_band_info(path)[0][b1-1][R]
    band_2 = get_band_info(path)[0][b2-1][R]
    kpoints = KPOINTS(path,fileformat='PROCAR')
    new_kpoints = kpoints.loc[R].reset_index()
    if ref_kpt == None:
        ref_kpt = new_kpoints.loc[0]
    k_norm = np.zeros(len(R))
    for k in ['kx','ky','kz']:
        k_norm += ((new_kpoints[k]-ref_kpt[k]))**2
    k_norm = np.sqrt(k_norm)
    if E0 == None :
        fit = FIT['GX']
        cond_E0 = False
    else :
        cond_E0 = True
        fit = lambda x,g,g2 : FIT['GX'](x,g,g2,E0)
    # GX :
    popt_GX_norm, err_GX_norm, a1_GX_norm, ea1_GX_norm, a2_GX_norm, ea2_GX_norm = rashba_p4mm(path, b1,b2,R_GX, fit, ref_kpt, *A, display = display)#ax = A[0], ay = A[1], az = A[2])
    popt_GX_ua, err_GX_ua, a1_GX_ua, ea1_GX_ua, a2_GX_ua, ea2_GX_ua = rashba_p4mm(path, b1,b2,R_GX, fit, ref_kpt, display = False)
    #popt_GX_norm, err_GX_norm, a1_GX_norm, ea1_GX_norm, a2_GX_norm, ea2_GX_norm = rashba_p4mm(path, b1,b2,R_GX, fit, ref_kpt, *A, display = display)#ax = A[0], ay = A[1], az = A[2])
    #popt_GX_ua, err_GX_ua, a1_GX_ua, ea1_GX_ua, a2_GX_ua, ea2_GX_ua = rashba_p4mm(path, b1,b2,R_GX, fit, ref_kpt, display = False)
    if display :
        print("              --------------------- GX ---------------------")
        print("fit GX : dE = 2*g*kx + 2*g2*kx**3 + E0")
    if cond_E0 :
        g, g2= popt_GX_norm
        eg, eg2 = err_GX_norm
        eE0 = 0
        g2_norm = g2
        eg2_norm = eg2
    else :
        g, g2, E0 = popt_GX_norm
        eg, eg2, eE0 = err_GX_norm
        g2_norm = g2
        eg2_norm = eg2
    if display :
        a1 =  a1_GX_norm
        ea1 = ea1_GX_norm
        a2 =  a2_GX_norm
        ea2 = ea2_GX_norm
        print("Normalisé : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g2 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2, g,eg,g2,eg2,E0,eE0))
    if cond_E0 :
        g, g2 = popt_GX_ua
        eg, eg2 = err_GX_ua
        eE0 = 0
        g2_ua = g2
        eg2_ua = eg2
    else :
        g, g2, E0 = popt_GX_ua
        eg, eg2, eE0 = err_GX_ua
        g2_ua = g2
        eg2_ua = eg2
        
    if display :
        a1 =  a1_GX_ua
        ea1 = ea1_GX_ua
        a2 =  a2_GX_ua
        ea2 = ea2_GX_ua
        print("ua        : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g2 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2,g,eg,g2,eg2,E0,eE0))
        print("              ----------------------------------------------")
    #plt.plot(k_norm, band_2-band_1, color = 'black')
    #plt.plot(k_norm, FIT['GX'](k_norm, *popt_GX_ua,E0),color = 'red', linestyle = ':')
    # GM :
    
    popt_GM_norm, err_GM_norm, a1_GM_norm, ea1_GM_norm, a2_GM_norm, ea2_GM_norm = rashba_p4mm(path, b1,b2,R_GM, fit, ref_kpt, *A, display = display)
    popt_GM_ua, err_GM_ua, a1_GM_ua, ea1_GM_ua, a2_GM_ua, ea2_GM_ua = rashba_p4mm(path, b1,b2,R_GM, fit, ref_kpt, display = False)
    
    if display :
        print("              --------------------- GM ---------------------")
        print("fit GM : dE = 2*g*kx + 2*(g2+g3)*kx**3 + E0")
    if cond_E0 :
        g, gt = popt_GM_norm
        err_GM_norm[1] += err_GX_norm[1]
        eg, egt = err_GM_norm
        eE0 = 0
    else :
        g, gt, E0 = popt_GM_norm
        err_GM_norm[1] += err_GX_norm[1]
        eg, egt, eE0 = err_GM_norm
    if display :
        a1 =  a1_GM_norm
        ea1 = ea1_GM_norm
        a2 =  a2_GM_norm
        ea2 = ea2_GM_norm
        print("Normalisé : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g1 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2,g,eg,gt-g2_norm,egt,E0,eE0))
    if cond_E0 :
        g, gt = popt_GM_ua
        err_GM_ua[1] += err_GX_ua[1]
        eg, egt = err_GM_ua
        eE0 = 0
    else :
        g, gt, E0 = popt_GM_ua
        err_GM_ua[1] += err_GX_ua[1]
        eg, egt, eE0 = err_GM_ua
    if display :
        a1 =  a1_GM_ua
        ea1 = ea1_GM_ua
        a2 =  a2_GM_ua
        ea2 = ea2_GM_ua
        print("ua        : alpha_1 = %.4f eV.A**2 (%.2f); alpha_2 = %.4f eV.A**2 (%.2f); g= %.4f eV.A (%.2f); g1 = %.4f eV.A**3 (%.2f); E0 = %.2f eV (%.2f)" %(a1,ea1,a2,ea2,g,eg,gt-g2_ua,egt,E0,eE0))
        print("              ----------------------------------------------")
    RES = {}
    RES['norm'] ={'GX' :  {'a1' : a1_GX_norm , 'a2' : a2_GX_norm, 'g': popt_GX_norm[0], 'g2' : popt_GX_norm[1]}, 
                  'GM' :  {'a1' : a1_GM_norm , 'a2' : a2_GM_norm, 'g': popt_GM_norm[0], 'g1' : popt_GM_norm[1]}}
    RES['ua'] ={'GX' :    {'a1' : a1_GX_ua   , 'a2' : a2_GX_ua,   'g': popt_GX_ua[0],   'g2' : popt_GX_ua[1]}, 
                'GM' :    {'a1' : a1_GM_ua   , 'a2' : a2_GM_ua,   'g': popt_GM_ua[0],   'g1' : popt_GM_ua[1]}}
    RES['error'] ={'GX' : {'a1' : ea1_GX_ua  , 'a2' : ea2_GX_ua,  'g': err_GX_ua[0],    'g2' : err_GX_ua[1]}, 
                   'GM':  {'a1' : ea1_GM_ua  , 'a2' : ea2_GM_ua,  'g': err_GM_ua[0],    'g1' : err_GM_ua[1]}}
    return RES
"""

def comp_fit(path, b1, b2, R_GX, R_GM, Ef = 0, E0 = None, ref_kpt = None, AX = None, args_aff = {}, display = True, opposite = True,
             g = None, g1 = None, g2 = None):
    RES = rashba_auto_p4mm(path, b1,b2,R_GX,R_GM, ref_kpt = ref_kpt, E0=E0, display = display, g = g, g1 = g1, g2 = g2)
    if AX == None :
        ax_GX = plt.subplot(2,2,1)
        ax_GX_dE = plt.subplot(2,2,2)
        ax_GM = plt.subplot(2,2,3)
        ax_GM_dE = plt.subplot(2,2,4)
    else :
        ax_GX = AX[0]
        ax_GX_dE = AX[1]
        ax_GM = AX[2]
        ax_GM_dE = AX[3]
    # --------- DIRECTION GX ---------
    g = RES['ua']['GX']['g']
    g2 = RES['ua']['GX']['g2']
    if E0 == None :
        E0 = 0
    fit = lambda x : FIT['GX'](x,g,g2,E0)
    R = R_GX
    kpoints = KPOINTS(path,fileformat='PROCAR')
    band_1 = get_band_info(path)[0][b1-1][R]-Ef
    band_2 = get_band_info(path)[0][b2-1][R]-Ef
    fit_band = lambda x,a : a*x**2 + g*x + g2*x**3 + band_1[0]
    
    new_kpoints = kpoints.loc[R].reset_index()
    if ref_kpt == None:
        ref_kpt = new_kpoints.loc[0]
    k_norm = np.zeros(len(R))
    for k in ['kx','ky','kz']:
        k_norm += ((new_kpoints[k]-ref_kpt[k]))**2
    k_norm = np.sqrt(k_norm)
    
    
    a1 = RES['ua']['GX']['a1']
    a2 = RES['ua']['GX']['a2']
    
    fit_band = lambda x,a : a*x**2 - (g*x + g2*x**3) + band_1[0]
    ax_GX.plot(k_norm, band_1, color ='red')
    ax_GX.plot(k_norm, band_2, color ='blue')
    ax_GX.plot(k_norm, fit_band(k_norm, a1), color ='red', linestyle = ':')
    fit_band = lambda x,a : a*x**2 + (g*x + g2*x**3) + band_2[0]
    ax_GX.plot(k_norm, fit_band(k_norm, a1), color ='blue', linestyle = ':')    
    ax_GX_dE.plot(k_norm, band_2-band_1, color ='black')
    ax_GX_dE.plot(k_norm, fit(k_norm), color ='green', linestyle = ':')
    
    if opposite == True :
        fit_band = lambda x,a : a*x**2 - (g*x + g2*x**3) + band_1[0]
        ax_GX.plot(-k_norm, band_1, color ='red')
        ax_GX.plot(-k_norm, band_2, color ='blue')
        ax_GX.plot(-k_norm, fit_band(-k_norm, a1), color ='red', linestyle = ':')
        fit_band = lambda x,a : a*x**2 + (g*x + g2*x**3) + band_2[0]
        ax_GX.plot(-k_norm, fit_band(-k_norm, a1), color ='blue', linestyle = ':')    
        ax_GX_dE.plot(-k_norm, band_2-band_1, color ='black')
        ax_GX_dE.plot(-k_norm, fit(k_norm), color ='green', linestyle = ':')
    
    # --------------------------------
    
    # --------- DIRECTION GM ---------
    
    g = RES['ua']['GM']['g']
    g1 = RES['ua']['GM']['g1']
    gt = g2 + g1
    fit = lambda x : FIT['GM'](x,g,g1+g2,E0)
    #fit = lambda x : FIT['GX'](x,g,gt,E0)
    R = R_GM
    kpoints = KPOINTS(path,fileformat='PROCAR')
    band_1 = get_band_info(path)[0][b1-1][R]-Ef
    band_2 = get_band_info(path)[0][b2-1][R]-Ef    
    new_kpoints = kpoints.loc[R].reset_index()
    if ref_kpt == None:
        ref_kpt = new_kpoints.loc[0]
    k_norm = np.zeros(len(R))
    for k in ['kx','ky','kz']:
        k_norm += ((new_kpoints[k]-ref_kpt[k]))**2
    k_norm = np.sqrt(k_norm)
    
    a1 = RES['ua']['GM']['a1']
    a2 = RES['ua']['GM']['a2']
    fit_band = lambda x,a : a*x**2 - (g*x + 1/2*(g1+g2)*x**3) + band_1[0]
    ax_GM.plot(k_norm, band_1, color ='red')
    ax_GM.plot(k_norm, band_2, color ='blue')
    ax_GM.plot(k_norm, fit_band(k_norm, a1), color ='red', linestyle = ':')
    fit_band = lambda x,a : a*x**2 + (g*x + 1/2*(g1+g2)*x**3) + band_2[0]
    ax_GM.plot(k_norm, fit_band(k_norm, a1), color ='blue', linestyle = ':')    
    ax_GM_dE.plot(k_norm, band_2-band_1, color ='black')
    ax_GM_dE.plot(k_norm, fit(k_norm), color ='green', linestyle = ':')
    
    if opposite == True :
        fit_band = lambda x,a : a*x**2 - (g*x + 1/2*(g1+g2)*x**3) + band_1[0]
        ax_GM.plot(-k_norm, band_1, color ='red')
        ax_GM.plot(-k_norm, band_2, color ='blue')
        ax_GM.plot(-k_norm, fit_band(-k_norm, a1), color ='red', linestyle = ':')
        fit_band = lambda x,a : a*x**2 + (g*x + 1/2*(g1+g2)*x**3) + band_2[0]
        ax_GM.plot(-k_norm, fit_band(-k_norm, a1), color ='blue', linestyle = ':')  
        ax_GM_dE.plot(-k_norm, band_2-band_1, color ='black')
        ax_GM_dE.plot(-k_norm, fit(k_norm), color ='green', linestyle = ':')
    
    
    
    #ax_GM.plot(-k_norm, band_2, color ='blue')
    #ax_GM.plot(-k_norm, fit_band(k_norm, *popt), color ='red', linestyle = ':')
    #ax_GM_dE.plot(-k_norm, band_2-band_1, color ='black')
    #ax_GM_dE.plot(-k_norm, fit(k_norm), color ='green', linestyle = ':')
    # --------------------------------
    
    return RES
    
    

def plot_ST_and_BS(path_BS, path_ST,b1,b2,R1 ,R2,Ef = None, E = None, 
                   E_dn = None, E_up = None,
                   g = None,g1 = None,g2 = None,
klim = [-0.5,0.5],

figname = 'ST',
fs_tick = 12,
fs_tick_bs = None,
tick_width = 1.5,
fs_suptitle = 14 ,
fs_title = 12,
fs_title_bs = None,
fs_contour = 10,
display_Spin_text = True, display_BS = True, display_info = True,
unit = 'norm',
quiver_space = 2,
quiver_space_tot = None,
quiver_head = 4,
quiver_width = 0.003,
quiver_scale = 20,
quiver_scale_tot = 30,
show_colorbar = True,
cmap = 'gist_rainbow',
show_contour = True,
HS = {'GX' : 'ZR', 'GM' : 'ZA'},
        #HS = {'GX' : 'GX', 'GM' : 'GM'}
ref_kpt = {'kx' : 0, 'ky' : 0, 'kz' : 0.5},
dE0 = 0):
    path = path_ST
    if quiver_space_tot == None:
        quiver_space_tot = quiver_space
    if Ef == None : 
        Ef = get_Ef_from_procar(path_BS)
    if fs_title_bs == None :
        fs_title_bs = fs_title
    if fs_tick_bs == None :
        fs_tick_bs = fs_tick
    cmap = plt.get_cmap(cmap)
    cmap.set_under('k',alpha=0)
    cmap.set_over('k',alpha=0)
    
    kpoints = KPOINTS(path_BS,fileformat='PROCAR')
    band_1 = get_band_info(path_BS)[0][b1-1]
    band_2 = get_band_info(path_BS)[0][b2-1]
    
    fig = plt.figure(figname,figsize = (16,6))
    #grid = plt.GridSpec(1, 2, wspace=0.01, hspace=0.4, left = 0.05,right = 0.95)
    #ax_ST = plt.subplot(grid[0])
    #plt.gca().set_aspect('equal', adjustable='box')
    #ax_ST_tot = plt.subplot(grid[0])
    #plt.gca().set_aspect('equal', adjustable='box')
    grid = plt.GridSpec(2, 8, wspace=0.5, hspace=0.4, left = 0.05,right = 0.95)
    if display_BS and display_Spin_text :
        grid = plt.GridSpec(2, 8, wspace=0.5, hspace=0.4, left = 0.05,right = 0.95)
    elif display_BS :
        grid = plt.GridSpec(2, 2, wspace=0.2, hspace=0.4, left = 0.05,right = 0.95)
    if display_BS :
        plt.rc('font', size=fs_tick_bs/1.5)
        ax_ZR = plt.subplot(grid[0,0])
        ax_ZR.set_title('ZR BS', fontsize = fs_title_bs)
        plt.rc('font', size=fs_tick_bs/1.5)
        
        ax_ZR_dE = plt.subplot(grid[0,1])
        ax_ZR_dE.set_title('ZR dE', fontsize = fs_title_bs)
        plt.rc('font', size=fs_tick_bs/1.5)
        
        ax_ZA = plt.subplot(grid[1,0])
        ax_ZA.set_title('ZA BS', fontsize = fs_title_bs)
        plt.rc('font', size=fs_tick_bs/1.5)
        
        ax_ZA_dE = plt.subplot(grid[1,1])
        ax_ZA_dE.set_title('ZA dE', fontsize = fs_title_bs)
        plt.rc('font', size=fs_tick_bs/1.5)
        for ax in [ax_ZR,ax_ZR_dE,ax_ZA,ax_ZA_dE]:
            ax.tick_params(labelsize = fs_tick_bs)
            ax.tick_params(labelsize = fs_tick_bs)
            ax.minorticks_on()
            ax.ticklabel_format(style='sci', axis='both',scilimits=(-2,3))
            ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick_bs, width = tick_width)
            ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick_bs/2)
            #ax.set_xlim(*klim)
            
        AX = [ax_ZR,ax_ZR_dE,ax_ZA,ax_ZA_dE]
        RES = comp_fit(path_BS, b1,b2,R1,R2, ref_kpt = ref_kpt, E0 = dE0, Ef = Ef, display = display_info, AX = AX, g = g, g1 = g1, g2 =g2)
        
        ax_ZR.set_ylim(*ax_ZR.get_ylim())
        ax_ZR.vlines(0,*ax_ZR.get_ylim(), linestyle = ':', color = 'black')
        
        ax_ZR_dE.set_ylim(*ax_ZR_dE.get_ylim())
        ax_ZR_dE.vlines(0,*ax_ZR_dE.get_ylim(), linestyle = ':', color = 'black')
        
        ax_ZA.set_ylim(*ax_ZA.get_ylim())
        ax_ZA.vlines(0,*ax_ZA.get_ylim(), linestyle = ':', color = 'black')
        
        ax_ZA_dE.set_ylim(*ax_ZA_dE.get_ylim())
        ax_ZA_dE.vlines(0,*ax_ZA_dE.get_ylim(), linestyle = ':', color = 'black')
        
        string = ''
        units = {'a1' : 'eV.A^2','a2' : 'eV.A^2','g' : 'eV.A','g1' : 'eV.A^3','g2' : 'eV.A^3'}
        for k in RES[unit].keys():
            string +="%2s :" %(k)
            for coef in RES[unit][k].keys():
                string += " %3s = %2.2e %6s (%.1f " %(coef,RES[unit][k][coef], units[coef], RES['error'][k][coef]) + '%)'
                if not coef == list(RES[unit][k].keys())[-1]:
                    string += ';'
            string +='\n'
        fig.suptitle(string, fontsize = fs_title)
        
        
        
    
    
    # -------------------------------------------------------------------- #
    # ------------------------ Texture de Spin --------------------------- #
    # -------------------------------------------------------------------- #
    if display_Spin_text :
        st = ST(path,b1)
        st.E -= Ef
        st2 = ST(path,b2)
        st2.E -= Ef
        if E_up == None :
            E_up = np.max(st.E)
        if E_dn == None :
            E_dn = np.min(st.E)
        Elim = E_up, E_dn
        e = np.linspace(E_dn, E_up,11)
        norm = plt.Normalize(E_dn, E_up)
        
        if E == None :
            E = ((E_up + E_dn )/2+E_dn)/2
        ax_ST = plt.subplot(grid[0:,2:5])
        plt.gca().set_aspect('equal', adjustable='box')
        ax_ST_tot = plt.subplot(grid[0:,5:8])
        plt.gca().set_aspect('equal', adjustable='box')
        for ax in [ax_ST,ax_ST_tot]:
            ax.tick_params(labelsize = fs_tick)
            ax.minorticks_on()
            ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick, width = tick_width)
            ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick/2)
            ax.set_xlim(*klim)
            ax.set_ylim(*klim)
        
        ax_ST.set_title('cut at %.2f' %(E), fontsize = fs_title)
        st.plot_ST_Ecut(E, color = 'red',quiver_space = quiver_space,
                        headlength=quiver_head,headwidth=quiver_head, width = quiver_width,scale = quiver_scale, ax = ax_ST)
        st2.plot_ST_Ecut(E, color = 'blue',quiver_space = quiver_space,
                         headlength=quiver_head,headwidth=quiver_head, width = quiver_width,scale = quiver_scale, ax = ax_ST)
        ax_ST.vlines(0,*ax_ST.get_ylim(), linestyle = ':', color = 'black')
        ax_ST.hlines(0,*ax_ST.get_xlim(), linestyle = ':', color = 'black')
        
        #ax_ST_tot = plt.subplot(grid[1])
        ax_ST_tot.set_title('total spin texture', fontsize = fs_title)
        
        
        
        
        l = st.plot_ST(Ef = 0, ax = ax_ST_tot,quiver_space = quiver_space_tot,
                       headlength=quiver_head,headwidth=quiver_head, width = quiver_width,scale = quiver_scale_tot, 
                       cmap = cmap, norm = norm)
        if show_contour :
            CS = ax_ST_tot.contour(st.kx,st.ky,st.E[:,:], e, cmap = cmap, linestyles = ':', linewidths = 1, norm = norm)
            ax_ST_tot.clabel(CS, inline=1, fontsize=10, fmt = '%.2f eV')
        ax_ST_tot.vlines(0,*ax_ST_tot.get_ylim(), linestyle = ':', color = 'black')
        ax_ST_tot.hlines(0,*ax_ST_tot.get_xlim(), linestyle = ':', color = 'black')
        
        bbox = ax_ST_tot.get_position().get_points()
        rect = [bbox[1][0], bbox[0][1], 0.007, bbox[1][1]-bbox[0][1]]
        if show_colorbar :
            cax = fig.add_axes(rect)
            cax.tick_params(labelsize = fs_tick)
            fig.colorbar(l,cax = cax)
            cax.set_ylabel('E-Ef (eV)', rotation = 270, va = 'bottom', fontsize = fs_title)
    print(quiver_space)
    if display_BS :
        return fig, RES
    else :
        return fig, None
#fig.tight_layout()
    HS = {'GX' : 'XG', 'GM' : 'XM'}
    ref_kpt = {'kx' : 0.5, 'ky' : 0.0, 'kz' : 0}
    columns = ['lattice', 'outplane','etaxx', 'g(XM)','g_error(XM)', 'g2(XM)', 'g2_error(XM)',
               'g(XG)', 'g_error(XG)','g1(XG)', 'g1_error(XG)']
    columns_2 = ['lattice',  'outplane','etaxx', 'g(MX)', 'g2(MX)', 
                 'g(MG)', 'g1(MG)']
    columns_2 = ['lattice',  'outplane','etaxx', 'g(XM)', 'g2(XM)', 
                 'g(XG)', 'g1(XG)']
    


def plot_ST_tot(path, band, sup = None, Ef = 0, Emin = 0, Emax = 0):
    quiver_space_tot = 1
    fs_title = 12
    quiver_head = 4
    quiver_width = 0.006
    quiver_scale = 20
    quiver_scale_tot = 60
    fs_tick = 12
    tick_width = 1.5
    show_colorbar = True
    st = ST(path,band)
    st.E -= Ef
    fig = plt.figure(sup, figsize = (4,3))
    grid = plt.GridSpec(1, 20, wspace=0.1, hspace=1, left = 0.2,right = 0.8, bottom = 0.2)
    ax = plt.subplot(grid[0:19])
    
    
    #ax.set_title(sup)
    ax.tick_params(labelsize = fs_tick)
    ax.minorticks_on()
    ax.ticklabel_format(style='sci', axis='both',scilimits=(-2,3))
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick, width = tick_width)
    ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick/2)
    ax.set_xlabel(r'$k_x$')
    ax.set_ylabel(r'$k_y$')
    plt.gca().set_aspect('equal', adjustable='box')
    norm = plt.Normalize(np.min(st.E),np.max(st.E))
    norm = plt.Normalize(Emin,Emax)
    
    cmap = plt.get_cmap('brg')
    l = st.plot_ST(ax = ax,quiver_space = quiver_space_tot,
                       headlength=quiver_head,headwidth=quiver_head, width = quiver_width,scale = quiver_scale_tot, 
                       cmap = cmap, norm=norm)
    ax.set_xlim(*ax.get_xlim())
    ax.set_ylim(*ax.get_ylim())
    klim = [-0.15,0.15]
    ax.set_xlim(*klim)
    ax.set_ylim(*klim)
    ax.vlines(0,*ax.get_ylim(), linestyle = ':', color = 'black')
    ax.hlines(0,*ax.get_xlim(), linestyle = ':', color = 'black')
    bbox = ax.get_position().get_points()
    rect = [bbox[1][0], bbox[0][1], 0.01, bbox[1][1]-bbox[0][1]]
    if show_colorbar :
        #cax = fig.add_axes(rect)
        cax = plt.subplot(grid[-1])
        cax.tick_params(labelsize = fs_tick)
        fig.colorbar(l,cax = cax)
        cax.set_ylabel(r'$E-\mathrm{E}_\mathrm{F}$'+' (eV)', rotation = 270, va = 'bottom', fontsize = fs_title)
    #cbar = plt.colorbar(l)
    #cbar.set_label('mz (mub)', rotation = 270)  

def plot_truc(path, band, sup = ''):
    
    quiver_space_tot = 3
    quiver_head = 6
    quiver_width = 0.003
    quiver_scale = 20
    quiver_scale_tot = 30
    fs_tick = 12
    tick_width = 1.5
    st = ST(path,band)
    fig = plt.figure(sup)
    ax = plt.subplot()
    
    
    ax.set_title(sup)
    ax.tick_params(labelsize = fs_tick)
    ax.minorticks_on()
    ax.ticklabel_format(style='sci', axis='both',scilimits=(-2,3))
    ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick, width = tick_width)
    ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick/2)
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    plt.gca().set_aspect('equal', adjustable='box')
    norm = plt.Normalize(-0.3,0.3)
    cmap = cmap = plt.get_cmap('gnuplot')
    l = st.plot_ST_contrib_z(ax = ax,quiver_space = quiver_space_tot,
                       headlength=quiver_head,headwidth=quiver_head, width = quiver_width,scale = quiver_scale_tot, 
                       cmap = cmap, norm = norm)
    cbar = plt.colorbar(l)
    cbar.set_label('mz (mub)', rotation = 270)
    

def plot_second_truc(path, E, Ef, sup = '', colors = {}, num_bz = [1,1], xlim = None, ylim = None, 
                     quiver_space = 5,
                    quiver_head = 4,
                    quiver_width = 0.010,
                    quiver_scale = 10,
                    quiver_scale_tot = 30,
                    fs_tick = 12,
                    tick_width = 1.5, name = None, ML = 0.2):
    if name == None :
        fig = plt.figure(sup, figsize = (4,3))
    else :
        fig = plt.figure(name, figsize = (4,3))
        
    ax = plt.subplot()
    if xlim == None :
        ax.set_xlim(-0.5,-0.5+num_bz[0])
    else :
        ax.set_xlim(*xlim)
    if ylim == None :
        ax.set_ylim(-0.5,-0.5+num_bz[1])
    else :
        ax.set_ylim(*ylim)
    ax.set_title(sup)
    ax.tick_params(labelsize = fs_tick)
    ax.minorticks_on()
    ax.ticklabel_format(style='sci', axis='both',scilimits=(-2,3))
    ax.yaxis.set_major_locator(MultipleLocator(ML))
    ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick, width = tick_width)
    ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick/2)
    ax.set_xlabel(r'$k_x$')
    ax.set_ylabel(r'$k_y$')
    plt.gca().set_aspect('equal', adjustable='box')
    for band in get_bands_limit(path,[E,E], Ef = Ef):
        st = ST(path,band)
        if band in colors.keys() :
            color = colors[band]
        else :
            color = None
        st.plot_ST_Ecut(E,Ef,ax = ax, num_bz = num_bz,
                        quiver_space = quiver_space, headwidth = quiver_head, width = quiver_width,
                        scale = quiver_scale, color = color, args_cs={'linestyles':'-'})
    for i in range(num_bz[1]*2):
        ax.hlines(i/2, *ax.get_xlim(), linestyle = ':')
    for j in range(num_bz[0]*2):
        ax.vlines(j/2, *ax.get_ylim(), linestyle = ':')
    fig.tight_layout()
    
def plot_trois_truc(path, kpt, knorm,Elim, Ef, colors = {}, sup = '',
                     quiver_space = 3,
                    quiver_head = 6,
                    quiver_width = 0.003,
                    quiver_scale = 20,
                    quiver_scale_tot = 30,
                    fs_tick = 12,
                    tick_width = 1.5,figsize = None,
                    **kwargs):
    fig = plt.figure(sup, figsize = figsize)
    ax = plt.subplot()
    ax.tick_params(labelsize = fs_tick)
    ax.minorticks_on()
    ax.ticklabel_format(style='sci', axis='both',scilimits=(-2,3))
    ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick, width = tick_width)
    ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick/2)
    ax.set_xlabel('k')
    ax.set_ylabel('E-EF(eV)')
    bands = get_bands_limit(path, Elim, Ef)
    for b in bands :
        if b in colors.keys() :
            color = colors[b]
        else :
            color = None
        bands = get_band_info(path)[0][b-1][kpt]
        plt.plot(knorm,bands-Ef, color = color, **kwargs)
    ax.set_ylim(*ax.get_ylim())
    ax.set_xlim(knorm[0],knorm[-1])
    ax.vlines(0, 4.2,4.8, linestyle = ':')
    ax.set_xticks([])
    fig.tight_layout()
    return fig




def plot_DFT_vs_model(path, b1,b2, GX,GM,path_ST, E_cut,ref_kpt = {'kx':0,'ky':0,'kz':0.5}, pos_E0 = 0,
                     quiver_space = 3,
                    quiver_head = 6,
                    quiver_width = 0.003,
                    quiver_scale = 20,
                    quiver_scale_tot = 30,
                    fs_tick = 12,
                    tick_width = 1.5,figsize = None,sup = None,color_band = None, color_model = 'green',
                    **kwargs):
    
    
    Ham = Rashba_model()
    
    kpoint = KPOINTS(path)
    RES = rashba_auto_p4mm(path,b1,b2,GX, GM, ref_kpt = ref_kpt, E0 = 0)
    KX_GX = kpoint['kx'][GX]
    KY_GX = kpoint['ky'][GX]
    knorm_GX = np.sqrt((KX_GX-ref_kpt['kx'])**2+(KY_GX-ref_kpt['ky'])**2)
    KX_GM = kpoint['kx'][GM]
    KY_GM = kpoint['ky'][GM]
    knorm_GM = np.sqrt((KX_GM-ref_kpt['kx'])**2+(KY_GM-ref_kpt['ky'])**2)
    band_1_GX = get_band_info(path)[0][b1-1][GX]
    band_1_GM = get_band_info(path)[0][b1-1][GM]
    band_2_GX = get_band_info(path)[0][b2-1][GX]
    band_2_GM = get_band_info(path)[0][b2-1][GM]
    
    a = (RES['ua']['GX']['a1']+RES['ua']['GX']['a2'])/2
    a1 = RES['ua']['GM']['a1']
    a2 = RES['ua']['GM']['a2']
    g=RES['ua']['GX']['g']
    g_gm=RES['ua']['GM']['g']
    g2=RES['ua']['GX']['g2']
    g1=RES['ua']['GM']['g1']
    
    
    function_GX = lambda k,a,g,g1,g2,signe=1 : a*k**2 + signe*(g*k+g2*k**3)
    function_GM = lambda k,a,g,g1,g2,signe=1 : a*k**2 + signe*(g*k+1/2*(g1+g2)*k**3)
    fig = plt.figure(sup, figsize = figsize)
    grid = plt.GridSpec(1, 5, wspace=1.2, hspace=1, left = 0.125,right = 0.95, bottom = 0.125)
    ax = plt.subplot(grid[0:2])
    ax.tick_params(labelsize = fs_tick)
    ax.minorticks_on()
    ax.ticklabel_format(style='sci', axis='both',scilimits=(-2,3))
    ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick, width = tick_width)
    ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick/2)
    ax.set_xlabel('k')
    ax.set_ylabel('E-EF(eV)')
    
    dknorm_GM = knorm_GM
    dknorm_GX = -knorm_GX
    
    ax.plot(-knorm_GX,band_1_GX-Ef, color = color_band)
    ax.plot(-knorm_GX,band_2_GX-Ef, color = color_band)
    
    ax.plot(dknorm_GM,band_1_GM-Ef, color = color_band)
    #plt.plot(knorm_GM,function_GM(knorm_GM,a1,g,g1,g2,signe = -1)+np.min(band_1_GM), color = 'red', linestyle = ':')
    ax.plot(dknorm_GM,band_2_GM-Ef, color = color_band)
    #plt.plot(knorm_GM,function_GM(knorm_GM,a1,g,g1,g2,signe = 1)+np.min(band_1_GM), color = 'blue', linestyle = ':')
    
    Ham.create_func(a1,g,g1,g2)
    print(band_1_GM[pos_E0]-Ef,np.array(Ham.Em(KX_GM,KY_GM))[pos_E0]+band_1_GM[pos_E0]-Ef)
    
    E0 = band_1_GM[pos_E0]
    
    ax.plot(dknorm_GM,Ham.Em(KX_GM-ref_kpt['kx'],KY_GM-ref_kpt['ky'])+E0-Ef, color = color_model, linestyle = ':')
    ax.plot(dknorm_GM,Ham.Ep(KX_GM-ref_kpt['kx'],KY_GM-ref_kpt['ky'])+E0-Ef, color = color_model, linestyle = ':')
    
    ax.plot(-knorm_GX,Ham.Em((KX_GX-ref_kpt['kx']),(KY_GX)-ref_kpt['ky'])+E0-Ef, color = color_model, linestyle = ':')
    ax.plot(-knorm_GX,Ham.Ep((KX_GX-ref_kpt['kx']),(KY_GX)-ref_kpt['ky'])+E0-Ef, color = color_model, linestyle = ':')
    xlim = [-np.max(knorm_GX),np.max(knorm_GM)]
    xlim = ax.get_xlim()
    ax.set_xlim(*xlim)
    ax.set_ylim(*ax.get_ylim())
    ax.vlines(0,*ax.get_ylim(), linestyle = '-')
    ax.hlines(E_cut,*ax.get_xlim(), linestyle = '--')
    ax1 = ax
    
    ax = plt.subplot(grid[2:])
    ax.tick_params(labelsize = fs_tick)
    ax.minorticks_on()
    ax.ticklabel_format(style='sci', axis='both',scilimits=(-2,3))
    ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = fs_tick, width = tick_width)
    ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = fs_tick/2)
    plt.gca().set_aspect('equal', adjustable='box')
    st1 = ST(path_ST, b1)
    st2 = ST(path_ST, b2)
    st1.plot_ST_Ecut(E_cut,Ef = Ef, ax = ax, color = color_band)
    st2.plot_ST_Ecut(E_cut,Ef = Ef, ax = ax, color = color_band)
    ax.set_xlabel('kx')
    ax.set_ylabel('ky')
    
    ax.set_xlim(*ax1.get_xlim())
    ax.set_ylim(*ax1.get_xlim())
    
    
    return fig
