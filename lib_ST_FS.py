#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 09:09:30 2019

@author: gosteau
"""

from lib_BS import *
from scipy.interpolate import griddata
import pickle as pkl
import bz2
import matplotlib.cm as cm
from matplotlib.colors import Normalize


class SpT():
    def __init__(self, path, band, Ef = None, E = 0, orbs_contrib = None, ST = False, interp = None, npts = 500):
        print( path)
        if Ef == None:
            Ef = get_Ef(path, same_dir = True)
        self.Ef = Ef
        self.E = E
        self.b = band
        self.interp = interp
        self.npts = npts
        self.ST = ST
        self.n_kpt = int(grep("NKPTS" , path + "OUTCAR")[0][0].split()[3])
        self.orbs = None
        if orbs_contrib != None :
            if type(orbs_contrib) == int :
                self.orbs = [str(orbs_contrib)+'_tot']
            if type(orbs_contrib) == str :
                self.orbs = [orbs_contrib]
            elif type(orbs_contrib) == list or type(orbs_contrib) == type(np.array([]))  :
                if type(orbs_contrib[0]) == int :
                    self.orbs = []
                    for at in orbs_contrib :
                        self.orbs += [str(at)+'_tot']
                elif type(orbs_contrib[0]) == str :
                    self.orbs = []
                    for orb in orbs_contrib :
                        self.orbs += [orb]
                        
            elif type(orbs_contrib) == dict :
                self.orbs = []
                for at in orbs_contrib['At']:
                    for orb in orbs_contrib['orb']:
                        self.orbs += [str(at)+'_'+orb]
                
        self.band = None
        self.path = path
        self.Kx = None
        self.Ky = None
        self.contrib = None
        self.mx = None
        self.band_sx = None
        self.band_sy = None
        self.my = None
        self.contour_plot = None
        self.get_k_cut()
        if self.orbs != None :
            self.get_contrib()
        if self.ST == True:
            self.get_spin()
    
    
        
    def get_k_cut(self):
        if type(self.band) == type(None) :
            HEADS = get_HEAD(path = self.path)
            BANDS_NAME = self.path + 'BANDS' + '/' + 'BAND_'+str(self.b)
            self.band = pd.read_csv(BANDS_NAME, delim_whitespace=True, names = HEADS,skiprows = 1,nrows = self.n_kpt)
            
        if self.interp != None :
            kxmin = self.band['kx'].min()
            kxmax = self.band['kx'].max()
            kymin = self.band['ky'].min()
            kymax = self.band['ky'].max()
            kx = np.linspace(kxmin,kxmax,self.npts)
            ky = np.linspace(kymin,kymax,self.npts)
            KX,KY = np.meshgrid(kx,ky)
            z = griddata((self.band['kx'], self.band['ky']), self.band['Energy']-self.Ef, (KX,KY), method = self.interp)
        else :
            #print(self.band[['kx','ky','Energy']].sort_values(['kx','ky']))
            hdf = self.band[['kx','ky','Energy']].sort_values(['kx','ky']).drop_duplicates().reset_index().pivot('kx','ky','Energy')
            X=hdf.columns.values
            Y=hdf.index.values
            z=hdf.values-Ef
            KX,KY = np.meshgrid(X, Y)
        fig = plt.figure('contour_tmp', figsize = (1,1))
        CS = plt.contour(KX, KY, z, [self.E], linewidths=0.5)
        paths = CS.collections[0].get_paths()
        verts = [xx.vertices for xx in paths]
        points = np.concatenate(verts)
        plt.close(fig)
        self.Kx = points[:,0]
        self.Ky = points[:,1]
        self.contour_plot = [KX,KY,z]
    
    def get_contrib(self):
        tmp = np.zeros(len(self.band))
        for orb in self.orbs :
            tmp += self.band[orb]
        if self.interp == None :
            method = 'linear'
        else :
            method = self.interp
        self.contrib = griddata((self.band['kx'], self.band['ky']), tmp, (self.Kx, self.Ky), method = method)
    
    def get_spin(self):
        if self.orbs == None:
            orbs = ['tot']
        else :
            orbs = self.orbs
        if type(self.band_sx) == type(None):
            HEADS = get_HEAD(path = self.path)
            BANDS_NAME = self.path + 'BANDS' + '/' + 'BAND_'+str(self.b)+'_sx'
            self.band_sx = pd.read_csv(BANDS_NAME, delim_whitespace=True, names = HEADS, skiprows = 1,nrows = self.n_kpt)
            BANDS_NAME = self.path + 'BANDS' + '/' + 'BAND_'+str(self.b)+'_sy'
            self.band_sy = pd.read_csv(BANDS_NAME, delim_whitespace=True, names = HEADS, skiprows = 1,nrows = self.n_kpt)
        
        #if orbs[0] != 'tot' :
        tmp_sx = np.zeros(len(self.band_sx))
        for orb in self.orbs :
            tmp_sx += self.band_sx[orb]  
        tmp_sy = np.zeros(len(self.band_sy))
        for orb in self.orbs :
            tmp_sy += self.band_sy[orb]
        #
        if self.interp == None :
            method = 'linear'
        else :
            method = self.interp
        self.mx = griddata((self.band_sx['kx'], self.band_sx['ky']), tmp_sx, (self.Kx, self.Ky), method = method)
        self.my = griddata((self.band_sy['kx'], self.band_sy['ky']), tmp_sy, (self.Kx, self.Ky), method = method)
    
    def update(self, orbs_contrib, E = 0):
        if orbs_contrib != None :
            if type(orbs_contrib) == int :
                self.orbs = [str(orbs_contrib)+'_tot']
            if type(orbs_contrib) == str :
                self.orbs = [orbs_contrib]
            elif type(orbs_contrib) == list or type(orbs_contrib) == type(np.array([]))  :
                if type(orbs_contrib[0]) == int :
                    self.orbs = []
                    for at in orbs_contrib :
                        self.orbs += [str(at)+'_tot']
                elif type(orbs_contrib[0]) == str :
                    self.orbs = []
                    for orb in orbs_contrib :
                        self.orbs += [orb]
                        
            elif type(orbs_contrib) == dict :
                self.orbs = []
                for at in orbs_contrib['At']:
                    for orb in orbs_contrib['orb']:
                        self.orbs += [str(at)+'_'+orb]
        if E != None:
            if self.E != E :
                self.E = E
                self.get_k_cut()
        self.get_contrib()
        if self.ST :
            self.get_spin()

    
    def plot(self, **kwargs):
        ax = kwargs.get('ax', None)
        norm = kwargs.get('norm', Normalize(0,1))
        colormap = kwargs.get('colormap', cm.seismic)
        color = kwargs.get('color', None)
        ST = kwargs.get('ST', True)
        contour_lw = kwargs.get('contour_lw', 1)
        contour_ls = kwargs.get('contour_ls', '-')
        quiver_head = kwargs.get('quiver_head', 4)
        quiver_width = kwargs.get('quiver_width', 0.003)
        quiver_scale = kwargs.get('quiver_scale', 10)
        quiver_space = kwargs.get('quiver_space', 1)
        inter_plot = kwargs.get('inter_plot', True)
        guide = kwargs.get('guide', True)
        if inter_plot == False :#and self.interp != None:
            KX,KY,z  = self.contour_plot
        else :
            kxmin = self.band['kx'].min()
            kxmax = self.band['kx'].max()
            kymin = self.band['ky'].min()
            kymax = self.band['ky'].max()
            kx = np.linspace(kxmin,kxmax,self.npts)
            ky = np.linspace(kymin,kymax,self.npts)
            KX,KY = np.meshgrid(kx,ky)
            z = griddata((self.band['kx'], self.band['ky']), self.band['Energy']-self.Ef, (KX,KY), method = 'linear')
            
        if color == None and self.orbs != None :
            
            color = colormap(norm(self.contrib))
        else :
            color = 'black'
        if ax == None :
            ax = plt.subplot()
        if self.orbs != None :
            contour_lw = 0
        if guide == True:
            CS = ax.contour(KX, KY, z, [self.E], colors = color, linewidths = contour_lw,linestyles = contour_ls)
        if self.orbs != None and guide == True :
            paths = CS.collections[0].get_paths()
            verts = [xx.vertices for xx in paths]
            points = np.concatenate(verts)
            KX = points[:,0]
            KY = points[:,1]
            ax.scatter(KX, KY, c = color, marker = '.', s = 0.1)
        if ST == True and self.ST == True:
            ax.quiver(self.Kx[::quiver_space],self.Ky[::quiver_space],self.mx[::quiver_space],self.my[::quiver_space], color = color,
                       scale_units='xy',angles = 'xy',headaxislength=quiver_head, width = quiver_width,scale = quiver_scale)
        


def plot_ST(path, used_band,Ef = None, E = 0,  figname= 'ST', ST = True, orb = None,interp = None, display_color = True,STs = {},  legend = False,**kwargs):
    fig = plt.figure(figname,figsize=(5,5))
    plt.title(figname)
    ax = plt.subplot()
    for b in used_band :
        if b in STs.keys() :
            st = STs[b]
            st.update(orb, E)
            st.plot(**kwargs)
        else :
            st = SpT(path,b,Ef,E, interp = interp, orbs_contrib = orb, ST = ST)
            st.plot(**kwargs)
            STs[b] = st
             
            
    plt.xlabel('kx')
    plt.ylabel('ky')
    
    xlim = kwargs.get('xlim', plt.xlim())
    ylim = kwargs.get('ylim', plt.ylim())
    plt.xlim(*xlim)
    plt.ylim(*ylim)
    plt.hlines(0,*xlim,linestyle=':', color = 'black')
    plt.vlines(0,*ylim,linestyle=':', color = 'black')
    
    norm = kwargs.get('norm', Normalize(0,1))
    colormap = kwargs.get('colormap', cm.seismic)
    color = kwargs.get('color', None)
    if color == None and orb != None and display_color and legend :
        cbar_ax = fig.add_axes([0.93, 0.1, 0.01, 0.8])
        mpl.colorbar.ColorbarBase(cbar_ax, cmap=colormap,
                                norm=norm)
    return STs, fig

#path = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/ISYM/"
Ef = 6.11877962
band = 43
E = -0.3
head = 4
width = 0.001
scale = 10

#st1 = SpT(path,band,Ef,E, interp = None, orbs_contrib = 'p')

ref = time.time()
#st1.get_k_cut()
#st1.get_contrib()
#st1.get_spin()
#print time.time()-ref

ref = time.time()
#st2.get_k_cut()
#st2.get_contrib()
#st2.get_spin()
#print time.time()-ref

colormap = cm.seismic
norm = Normalize(0,1)
#plot_ST(path,[44,45,46,47,48,49],Ef, 2.5, interp = 'cubic',quiver_scale = 15, orb = 'd', color = 'black')
#plt.scatter(st2.Kx,st2.Ky, color = colormap(norm(st2.contrib)), marker = '+')
#st2.get_k_cut(cache_fig = False, colors = colormap(norm(st2.contrib)))
"""
fig = plt.figure('ST',figsize=(10,10))
st2.plot()
plt.xlabel('kx')
plt.ylabel('ky')
xlim = plt.xlim()
plt.xlim(*xlim)
ylim = plt.ylim()
plt.ylim(*ylim)
plt.hlines(0,*xlim,linestyle=':', color = 'black')
plt.vlines(0,*ylim,linestyle=':', color = 'black')
"""
#plt.quiver(st2.Kx,st2.Ky,st2.mx,st2.my, color = colormap(norm(st2.contrib)),scale_units='xy',angles = 'xy',headaxislength=head, width = width,scale = scale)
C = np.array([[255,0,0]])
cmap_r = mpl.colors.ListedColormap(C/255.0)
cmap_r.set_under('k',alpha=0)
C = np.array([[0,255,0]])
cmap_g = mpl.colors.ListedColormap(C/255.0)
cmap_g.set_under('k',alpha=0)
C = np.array([[0,0,255]])
cmap_b = mpl.colors.ListedColormap(C/255.0)
cmap_b.set_under('k',alpha=0)
C = np.array([[0,0,0]])
cmap_n = mpl.colors.ListedColormap(C/255.0)
cmap_n.set_under('k',alpha=0)
C = np.array([[255,0,255]])
cmap_p = mpl.colors.ListedColormap(C/255.0)
cmap_p.set_under('k',alpha=0)
C = np.array([[0,255,255]])
cmap_c = mpl.colors.ListedColormap(C/255.0)
cmap_c.set_under('k',alpha=0)
b = displayed_bands(path,1.5,3.5,Ef)
#path = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/kz0/"
#path = None
#plt.scatter(st1.Kx,st1.Ky, c = colormap(norm(st1.contrib)), marker = '+')
#plt.quiver(st1.Kx,st1.Ky,st1.mx,st1.my, color = colormap(norm(st1.contrib)),scale_units='xy',angles = 'xy',headaxislength=head, width = width,scale = scale)

def get_b_test(path,Ef):
    HEADS = get_HEAD(path = path)
    n_kpt = int(grep("NKPTS" , path + "OUTCAR")[0][0].split()[3])
    n_bands = int(grep("NKPTS" , path + "OUTCAR")[0][0].split()[-1])
    log = pd.DataFrame(index = range(1,n_bands+1), columns = ['Emin','Emax'])
    for b in range(1,n_bands+1):
        BANDS_NAME = path + 'BANDS' + '/' + 'BAND_'+str(b)
        band = pd.read_csv(BANDS_NAME, delim_whitespace=True, names = HEADS, 
                          usecols = ['Energy'],skiprows = 1,nrows = n_kpt)
        #print b, band['Energy'].min()-Ef, band['Energy'].max()-Ef
        log.loc[b]['Emin'] = band['Energy'].min()-Ef
        log.loc[b]['Emax'] = band['Energy'].max()-Ef
    return log

def get_bands_set(path,Ef,E, **kwargs):
    log = kwargs.get('log', None)
    if type(log) == type(None):
        log = get_b_test(path,Ef)
    bands = list(log.loc[log['Emin'] < E].loc[log['Emax']> E].index)
    return bands


def routine_BS_d(path, Ef = None,at = 'tot', path_Ef = None, HS_dir = ['X','G','X'], interp_kind = 'cubic' , ax = None, Elim = [-0.5,0.5], seuil = 0.2):
    
    norm = Normalize(seuil-0.1,seuil)
    
    routine(path, Ef, orb = 'tot', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_n, norm = norm, ax = ax)
    routine(path, Ef, orb = 'dxy', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_r, norm = norm, ax = ax)
    routine(path, Ef, orb = 'dyz', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_g, norm = norm, ax = ax)
    routine(path, Ef, orb = 'dxz', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_p, norm = norm, ax = ax)
    routine(path, Ef, orb = 'dz2', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_b, norm = norm, ax = ax)
    routine(path, Ef, orb = 'x2-y2', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_c, norm = norm,ax = ax)
    
def routine_BS_d10(path, Ef = None,at = 'tot', path_Ef = None, HS_dir = ['X','G','X'], interp_kind = 'cubic' , ax = None, Elim = [-0.5,0.5], seuil = 0.2):
    
    norm = Normalize(seuil-0.1,seuil)
    
    routine(path, Ef, orb = 'tot', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_n, norm = norm, ax = ax)
    routine(path, Ef, orb = 'd', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_r, norm = norm, ax = ax)
    routine(path, Ef, orb = 'p', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_b, norm = norm, ax = ax)
    routine(path, Ef, orb = 's', HS_dir = HS_dir, Elim = Elim, interp_kind=interp_kind, cmap = cmap_g, norm = norm, ax = ax)
    
def routine_d(path, Ef, E, seuil = 0.2, dE = 1,figname = 'TS',used_band = None, STs = {}, **kwargs):
    #plot_ST(path,b[:-2],Ef, E, interp = None,quiver_scale = 15, orb = None, 
    #        figname = figname, display_color=False, color = 'black', ST = False, **kwargs)
    if used_band == None :
        b = get_bands_set(path,Ef,E, **kwargs)
    else :
        b = used_band
    print('dx, dyz : green ;',end='')
    plot_ST(path,b[:],Ef, E,quiver_scale = 15, orb = ['dxz','dyz'], STs = STs,
            figname = figname, display_color=False, colormap = cmap_g, norm = Normalize(seuil-0.1,seuil), **kwargs)
    print( 'dxy : red ;',end='')
    plot_ST(path,b[:],Ef, E,quiver_scale = 15, orb = ['dxy'], STs = STs,
            figname = figname, display_color=False, colormap = cmap_r, norm = Normalize(seuil-0.1,seuil), **kwargs)
    print( 'p, s : black ;',end='')
    plot_ST(path,b[:],Ef, E,quiver_scale = 15, orb = ['px','py','pz','s'], STs = STs,
            figname = figname, display_color=False, colormap = cmap_n, norm = Normalize(seuil-0.1,seuil), **kwargs)
    print('eg : blue ;',end='')
    sts , fig= plot_ST(path,b[:],Ef, E,quiver_scale = 15, orb = ['dz2','x2-y2'], STs = STs,
            figname = figname, display_color=False, colormap = cmap_b, norm = Normalize(seuil-0.1,seuil), **kwargs)
    #print 's : purple ;'
    #plot_ST(path,b[:],Ef, E, interp = None,quiver_scale = 15, orb = ['s'], STs = STs,
    #        figname = figname, display_color=False, colormap = cmap_p, norm = Normalize(seuil-0.1,seuil), **kwargs)
    return fig
    #plot_ST(path,b[:],Ef, E, interp = None,quiver_scale = 15, orb = ['px','py','pz'], 
    #        figname = figname, display_color=False, colormap = cmap_g, norm = Normalize(seuil-0.1,seuil), **kwargs)
    
def routine_10(path, Ef, E, seuil = 0.2, dE = 1,figname = 'TS',used_band = None, STs = {}, **kwargs):
    #plot_ST(path,b[:-2],Ef, E, interp = None,quiver_scale = 15, orb = None, 
    #        figname = figname, display_color=False, color = 'black', ST = False, **kwargs)
    print(path)
    if used_band == None :
        b = get_bands_set(path,Ef,E, **kwargs)
    else :
        b = used_band
    print('d red',end='')
    plot_ST(path,b[:],Ef, E, interp = None,quiver_scale = 15, orb = ['d'],  STs = STs,
            figname = figname, display_color=False, colormap = cmap_g, norm = Normalize(seuil-0.1,seuil), **kwargs)
    print( 'p : green')
    sts, fig = plot_ST(path,b[:],Ef, E, interp = None,quiver_scale = 15, orb = ['p'],  STs = STs,
            figname = figname, display_color=False, colormap = cmap_n, norm = Normalize(seuil-0.1,seuil), **kwargs)
    #print 'p : m'
    return fig
    #plot_ST(path,b[:],Ef, E, interp = None,quiver_scale = 15, orb = ['px','py','pz'], 
    #        figname = figname, display_color=False, colormap = cmap_g, norm = Normalize(seuil-0.1,seuil), **kwargs)
path = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/kz05/LORBIT_10/"
path = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/kz0/"
path_05 = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/kz05/"
#path_05 = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/kz05/zoom/"
path_BS = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/BS/"
path_RZA = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/BS/RZA/"
path_XGM = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/BS/XGM/"
log = get_b_test(path,Ef)
log_05 = get_b_test(path_05,Ef)

"""
path_GZ = '/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/GZ/kz0/'
#log_GZ = get_b_test(path_GZ,Ef)

path_GZ_05 = '/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/GZ/kz05/'
#log_GZ_05 = get_b_test(path_GZ_05,Ef)
"""
#path_Ef = "/home/gosteau/Documents/These/1ere année/LAOSTO/asym/PBESOL/8_STO/relax_higher/next/more_highter/SOC SCF/More_points/"
#path = "/home/gosteau/Documents/These/1ere année/LAOSTO/asym/PBESOL/8_STO/relax_higher/next/more_highter/SOC SCF/More_points/NSCF_ST/"
#path_BS = "/home/gosteau/Documents/These/1ere année/LAOSTO/asym/PBESOL/8_STO/relax_higher/next/more_highter/SOC SCF/More_points/NSCF_BS/"
#Ef = get_Ef(path_Ef, same_dir = True)
#log = get_b_test(path,Ef)
"""
STs = {}
for e in np.linspace(-0.09,0,2):
    routine_d(path, Ef, e, STs = STs, log = log, xlim = [-0.2,0.2], ylim = [-0.2,0.2], figname = 'TS_'+str(e), interp = 'cubic')
"""
STs = {}
figname = 'ST'
de = 0.1
b = get_bands_set(path,Ef,1.7, log = log_05)[0:-1]
custom_lines = []
qs = 30

e = 1.65
c = cmap_b
plot_ST(path_05,b[:],Ef, e, figname = figname,interp = None,quiver_scale = qs, orb = ['tot'],  STs = STs,
            display_color=False, colormap = c, norm = Normalize(0.5-0.1,0.5))
custom_lines += [Line2D([0], [0], color=c(0.), label = "%1.2f" %(e))]

e += de
c = cmap_r
plot_ST(path_05,b[:],Ef, e, figname = figname, interp = None,quiver_scale = qs, orb = ['tot'],  STs = STs,
            display_color=False, colormap = c, norm = Normalize(0.5-0.1,0.5))
custom_lines += [Line2D([0], [0], color=c(0.), label = "%1.2f" %(e))]

e += de
c = cmap_n
plot_ST(path_05,b[:],Ef, e, figname = figname, interp = None,quiver_scale = qs, orb = ['tot'],  STs = STs,
            display_color=False, colormap = c, norm = Normalize(0.5-0.1,0.5))
custom_lines += [Line2D([0], [0], color=c(0.), label = "%1.2f" %(e))]

e += de
c = cmap_g
plot_ST(path_05,b[:],Ef, e, figname = figname, interp = None,quiver_scale = qs, orb = ['tot'],  STs = STs,
            display_color=False, colormap = c, norm = Normalize(0.5-0.1,0.5))
custom_lines += [Line2D([0], [0], color=c(0.), label = "%1.2f" %(e))]

e += de
c = cmap_p
plot_ST(path_05,b[:],Ef, e, figname = figname, interp = None,quiver_scale = qs, orb = ['tot'],  STs = STs,
            display_color=False, colormap = c, norm = Normalize(0.5-0.1,0.5))
custom_lines += [Line2D([0], [0], color=c(0.), label = "%1.2f" %(e))]

plt.xlim(-0.2,0.2)
plt.ylim(-0.2,0.2)
plt.legend(handles = custom_lines)
