#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:43:59 2018

@author: gosteau
"""
# System
import copy
import os as os
import time
import sys
import shutil

# Matplotlib
import matplotlib as mpl
mpl.use('nbagg')
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import pickle as pkl
import bz2
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator, IndexLocator )
from matplotlib.colors import Normalize
from scipy.interpolate import griddata
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['ps.fonttype'] = 42

# Set the default fontsize
font = {'size'   : 22}
font = {'size'   : 12}
mpl.rc('font', **font)
#mpl.rcParams['xtick.labelsize'] = mpl.rcParams['font.size']*3/4
mpl.rcParams['ytick.labelsize'] = mpl.rcParams['xtick.labelsize']
mpl.rcParams['legend.fontsize'] = mpl.rcParams['xtick.labelsize']
#mpl.rc('text', usetex=True)

COLORS = {'Sr' : 'green', 'Ti' : 'blue', 'La' : 'purple', 'Al' : 'orange', 'O' : 'red', 'Fe' : 'orange', 'Ir' : 'purple'}



# Calculation
import numpy as np
import pandas as pd
from math import ceil
from pylab import polyfit
from scipy.optimize import curve_fit
import scipy.integrate as scp
import scipy.interpolate as inter
from cffi import FFI


def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r"):
    """
    Call in a loop to create terminal progress bar
    @params:
        iteration   - Required  : current iteration (Int)
        total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """
    percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
    filledLength = int(length * iteration // total)
    bar = fill * filledLength + '-' * (length - filledLength)
    print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
    # Print New Line on Complete
    if iteration == total: 
        print()

def search(tab, val):
    values = []
    for i in range(len(tab)) :
        if tab[i] == val :
            values.append(i)
    return values

def clear_path(path,space = False, term=True):
    if term == True :
        if path[-1] != '/' :
            new_path = path + '/'
        else :
            new_path = path
    else :
        new_path = path
    if ' ' in path and space == True :
        new_path = ''
        for l in path :
            if ' ' == l :
                new_path += '\ '
            else :
                new_path += l
    
    return new_path

def grep(string, FILENAME) :
    FILE=open(FILENAME,'r')
    res=[[],[]]
    count = 0
    if type(string) == str :
        for l in FILE :
            count += 1
            if string in l :
                res[0].append(l.replace("\n",""))
                res[1].append(count)
    elif type(string) == list:
        for l in FILE :
            count += 1
            for s in string :
                if s in l :
                    Bool = True
                else :
                    Bool = False
                    break
            if Bool == True:
                res[0].append(l.replace("\n",""))
                res[1].append(count)
        
    FILE.close()
    return res

def get_row(n, FILENAME,a=0) :
    FILE=open(FILENAME,'r')
    res=''
    count = 0
    
    while count != n :
        count += 1
        res = FILE.readline().replace("\n","")
        if count == n :
            for i in range(a):
                res += '\n' + FILE.readline().replace("\n","")
            break
    FILE.close()
    return res

def get_lsorbit(file):
    try : 
        lsorbit = grep('LSORBIT', file)[0][0].split()[2]
        if lsorbit == 'T':
            lsorbit = True
        else :
            lsorbit = False
    except :
        lsorbit = False
    return lsorbit

def get_maxis(file) :
    n_SC, n_CS = grep('SAXIS', file)[1]
    M_SC = []
    for n in [n_SC+2,n_SC+3,n_SC+4] :
        M_SC.append(np.float64(np.array(get_row(n,file).split())[::2]))
    M_CS = []
    for n in [n_CS+2,n_CS+3,n_CS+4] :
        M_CS.append(np.float64(np.array(get_row(n,file).split())[::2]))
    return np.array(M_SC), np.array(M_CS)

class INFO_PROCAR():
    def __init__(self,path) :
        f = open(path+'PROCAR','r')
        for line in f :
            if line[0] == '#' :
                #self.n_kpts = int(line.split()[3])
                tmp = line.split('#')
                tmp_kpt = tmp[1].split(':')
                tmp_bands = tmp[2].split(':')
                tmp_ions = tmp[3].split(':')
                
                self.n_bands = int(tmp_bands[-1].split()[0])
                self.n_ions = int(tmp_ions[-1].split()[0])
                #self.n_bands = int(line.split()[7])
                #self.n_ions = int(line.split()[-1])
            if line[0:3] == 'ion':
                self.orbs = np.array(line.split()[1:])
                break
        f.close()
        self.n_kpts = wc(path+'BANDS/kpts_tmp')

class INFO_OUTCAR():
    def __init__(self, path):
        self.path = path+'OUTCAR'
        self.get_lattice()
        
    def get_lattice(self):
        try :
            line = grep('direct lattice vectors', self.path)[1][-1]
            Matrix = np.float64(get_row(line+1,self.path,a=3).split()).reshape((3,6))
            self.direct_lattice = Matrix[:,[0,1,2]]
            self.rec_lattice = Matrix[:,[3,4,5]]
            
        except :
            print('error: cannot find lattice, trying other method')
            try :
                line = grep('A1', self.path)[1][0]
                tmp = get_row(line+0,self.path,a=3).replace('A1','').replace('A2','').replace('A3','').replace('=','').replace('(','').replace(',','').replace(')','')
                Matrix = np.float64(tmp.split()).reshape((3,3))
                self.direct_lattice = Matrix
                a,b,c = Matrix
                V = np.dot(np.cross(a,b),c)
                u = np.cross(b,c)/V
                v = np.cross(c,a)/V
                w = np.cross(a,b)/V
                self.rec_lattice = np.array([u,v,w])
            except :
                print('error: cannot find lattice')
                self.direct_lattice = None
                self.rec_lattice = None
        
        try :   
            line = grep('length of vectors', self.path)[1][-1]
            Matrix = np.float64(get_row(line+1,self.path,a=0).split())
            self.length_direct_vectors = Matrix[[0,1,2]]
            self.length_rec_vectors = Matrix[[3,4,5]]
        except :
            print('error: cannot length of lattice, trying other method')
            try :
                a,b,c = self.direct_lattice
                self.length_direct_vectors = np.array([np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c) ])
                u,v,w = self.rec_lattice
                self.length_rec_vectors = np.array([np.linalg.norm(u), np.linalg.norm(v), np.linalg.norm(w) ])
            
            except :
                print('error: cannot length of lattice')
                self.length_direct_vectors = None
                self.length_rec_vectors = None
            




def wc(file):
    f = open(file,'r')
    l= 0;
    for line in f:
        l += 1
    f.close()
    return l

def get_band_info(path, **kwargs):
    print('hello')
    spin = kwargs.get('spin', '')
    info = kwargs.get('info', INFO_PROCAR(path))
    if os.path.isfile(path+"BANDS/band_info"+spin+".csv") == True :
        #bands, occ = np.load(open(path+"BANDS/band_info.csv",'rb'))
        bands, occ = np.load(open(path+"BANDS/band_info"+spin+".csv",'rb'))
    else :
        if os.path.isfile(path+"BANDS/band_info_tmp"+spin) == True :
            if 'x' in spin or 'y' in spin or 'z' in spin :
                spin = ''
            bands = []
            occ = []
            #print("band_info_tmp"+spin)
            try :
                f = open(path+"BANDS/band_info_tmp"+spin,'r')
            except :
                f = open(path+"BANDS/band_info_tmp"+'_up','r')

            l = 0
            for line in f :
                sp = line.split()
                if l%(info.n_bands) == 0 :
                    bands.append([])
                    occ.append([])
                bands[int(l/info.n_bands)].append(float(sp[0]))
                occ[int(l/info.n_bands)].append(float(sp[1]))
                l += 1
            f.close()
            print(len(bands), l, info.n_bands)
            bands[int(l/info.n_bands)] = np.array(bands[int(l/info.n_bands)])
            occ[int(l/info.n_bands)] = np.array(occ[int(l/info.n_bands)])
            bands = np.transpose(np.array(bands))
            occ = np.transpose(np.array(occ))
            np.save(open(path+"BANDS/band_info"+spin+".csv",'wb'),np.array([bands,occ]))
        else :
            raise Exception("file %s does not exist." %(path+"BANDS/band_info"+spin))
    return bands, occ

         
def get_band_info(path, **kwargs):
    spin = kwargs.get('spin', '')
    info = kwargs.get('info', INFO_PROCAR(path))
    if os.path.isfile(path+"BANDS/band_info"+spin+".csv") == True :
        #bands, occ = np.load(open(path+"BANDS/band_info.csv",'rb'))
        bands, occ = np.load(open(path+"BANDS/band_info"+spin+".csv",'rb'))
    else :
        if 'x' in spin or 'y' in spin or 'z' in spin :
            spin = ''
        bands = []
        occ = []
        #print("band_info_tmp"+spin)
        try :
            f = open(path+"BANDS/band_info_tmp"+spin,'r')
        except :
            f = open(path+"BANDS/band_info_tmp"+'_up','r')
        
        l = 0
        for line in f :
            sp = line.split()
            if l%(info.n_bands) == 0 :
                bands.append([])
                occ.append([])
            bands[int(l/info.n_bands)].append(float(sp[0]))
            occ[int(l/info.n_bands)].append(float(sp[1]))
            l += 1
        f.close()
        bands = np.transpose(np.array(bands))
        occ = np.transpose(np.array(occ))
        np.save(open(path+"BANDS/band_info"+spin+".csv",'wb'),np.array([bands,occ]))
    return bands, occ

def get_band_info(path, fileformat = 'PROCAR', name = 'OUTCAR', **kwargs):
    if fileformat == 'PROCAR' :
        spin = kwargs.get('spin', '')
        info = kwargs.get('info', INFO_PROCAR(path))
        if os.path.isfile(path+"BANDS/band_info"+spin+".csv") == True :
            #bands, occ = np.load(open(path+"BANDS/band_info.csv",'rb'))
            bands, occ = np.load(open(path+"BANDS/band_info"+spin+".csv",'rb'))
        else :
            if 'x' in spin or 'y' in spin or 'z' in spin :
                spin = ''
            bands = []
            occ = []
            #print("band_info_tmp"+spin)
            try :
                f = open(path+"BANDS/band_info_tmp"+spin,'r')
            except :
                f = open(path+"BANDS/band_info_tmp"+'_up','r')

            l = 0
            for line in f :
                sp = line.split()
                if l%(info.n_bands) == 0 :
                    bands.append([])
                    occ.append([])
                bands[int(l/info.n_bands)].append(float(sp[0]))
                occ[int(l/info.n_bands)].append(float(sp[1]))
                l += 1
            f.close()
            bands = np.transpose(np.array(bands))
            occ = np.transpose(np.array(occ))
            np.save(open(path+"BANDS/band_info"+spin+".csv",'wb'),np.array([bands,occ]))
        return bands, occ
    else :
        lines = grep("band No.", clear_path(path)+name)[-1]
        nbands = int(grep("NBANDS", clear_path(path)+name)[0][0].split()[-1])
        nkpts = int(grep("NKPTS", clear_path(path)+name)[0][0].split()[3])
        print(nbands, nkpts)
        bands = np.zeros((nbands, nkpts))
        occ = np.zeros((nbands, nkpts))
        for k in range(nkpts) :
            
            tmp = np.float64(get_row(lines[k]+1, clear_path(path)+name, nbands).split()).reshape((nbands, 3))
            bands[:,k] = tmp[:,1]
            occ[:,k] = tmp[:,2]
        return bands, occ


def get_Ef_from_procar(path, spin = ''):
    path = clear_path(path)
    if os.path.isfile(path+"BANDS/band_info_tmp"+spin) == True :
        bands, occ = get_band_info(path, spin = spin)
        ix = np.where(occ != 0)[0]
        E = np.max(bands[ix])
        return E
    elif os.path.isfile(path+"BANDS/band_info_tmp_up") == True and os.path.isfile(path+"BANDS/band_info_tmp_dn") == True:
            E1 = get_Ef_from_procar(path,spin='_up')
            E2 = get_Ef_from_procar(path,spin='_dn')
            if E2>E1 :
                return E2
            else :
                return E1
    """
    Res = [pd.DataFrame({'occ' : occ[b],'E': bands[b]}) for b in range(len(occ))]
    E = 0
    for b in range(len(occ)):
        try :
            E = np.array(Res[b].loc[Res[b]['occ'] != 0].sort_values(['occ','E'], ascending=[True,False])['E'])[0]
        except :
            break
    return E
    """

def get_VBM(path, spin = '',fileformat = 'OUTCAR'):
    if fileformat == 'PROCAR' :
        try :
            bands, occ = get_band_info(path, spin = spin)
        except :
            E1 = get_VBM(path,spin='_up', fileformat=fileformat)
            E2 = get_VBM(path,spin='_dn', fileformat=fileformat)
            if E2>E1 :
                return E2
            else :
                return E1   
        Res = [pd.DataFrame({'occ' : occ[b],'E': bands[b]}) for b in range(len(occ))]
        list_e = Res[0].loc[Res[0]['occ'] == 0]['E']
        if len(list_e) == 0 :
            E = -9e99
        else :
            E = list_e.max()
        for b in range(1,len(occ)):
            list_e = Res[b].loc[Res[b]['occ'] != 0]['E']
            if len(list_e) == 0 :
                Eb = -9e99
            else :
                Eb = list_e.max()
            if Eb > E :
                E = Eb
        return E
    else :
        return get_VBM_OUTCAR(path)
    
def get_CBM(path, spin = '',fileformat = 'OUTCAR'):
    if fileformat == 'PROCAR' :
        try :
            bands, occ = get_band_info(path, spin = spin)
        except :
            E1,k1 = get_CBM(path,spin='_up', fileformat=fileformat)
            E2,k2 = get_CBM(path,spin='_dn', fileformat=fileformat)
            if E2>E1 :
                return E2,k2
            else :
                return E1,k1   
        Res = [pd.DataFrame({'occ' : occ[b],'E': bands[b]}) for b in range(len(occ))]
        list_e = Res[0].loc[Res[0]['occ'] == 0]['E']
        
        if len(list_e) == 0 :
            E = 9e99
        else :
            E = list_e.min()
        
        
        for b in range(1,len(occ)):
            list_e = Res[b].loc[Res[b]['occ'] == 0]['E']
            if len(list_e) == 0 :
                Eb = 9e99
            else :
                Eb = list_e.min()
            if Eb < E :
                E = Eb
            k = np.where(bands == E)[0][0]
            return E,k
    else :
        return get_CBM_OUTCAR(path)
    
def get_CBM(path, spin = '',fileformat = 'PROCAR'):
    if fileformat == 'PROCAR' :
        try :
            bands, occ = get_band_info(path, spin = spin)
        except :
            E1,k1 = get_CBM(path,spin='_up', fileformat=fileformat)
            E2,k2 = get_CBM(path,spin='_dn', fileformat=fileformat)
            if E2>E1 :
                return E2,k2
            else :
                return E1,k1   
        Res = [pd.DataFrame({'occ' : occ[b],'E': bands[b]}) for b in range(len(occ))]
        list_e = Res[0].loc[Res[0]['occ'] == 0]['E']
        
        if len(list_e) == 0 :
            E = 9e99
        else :
            E = list_e.min()
        
        
        for b in range(1,len(occ)):
            list_e = Res[b].loc[Res[b]['occ'] == 0]['E']
            if len(list_e) == 0 :
                Eb = 9e99
            else :
                Eb = list_e.min()
            if Eb < E :
                E = Eb
            k = np.where(bands == E)[0][0]
            return E,k
    else :
        return get_CBM_OUTCAR(path)
  

      

def get_Energy_from_OUTCAR(path, filename = 'OUTCAR'):
    cpath = clear_path(path)
    nline = grep('k-point     1', cpath+filename)[1][0]
    VALUES = []
    f = open(cpath+filename, 'r')
    i = 0
    for line in f :
        i+=1
        if i < nline :
            continue
        else :
            sp = line.split()
            if len(sp) != 3 :
                if len(sp) == 1:
                    break
                else :
                    continue
            else :
                if sp[0] != 'spin' :
                    VALUES.append(sp)
            
    f.close()
    return np.float64(VALUES)

def get_VBM_OUTCAR(path):
    VALUES = get_Energy_from_OUTCAR(path)
    return VALUES[VALUES[:,2] != 0,1].max()

def get_CBM_OUTCAR(path):
    VALUES = get_Energy_from_OUTCAR(path)
    return VALUES[VALUES[:,2] == 0,1].min()


def get_gap(path, spin = '',fileformat = 'OUTCAR') :
    VBM = get_VBM(path, spin,fileformat)
    CBM = get_CBM(path, spin,fileformat)
    return CBM-VBM

def get_pola(path):
    try : 
        return np.float64(grep('p[elc]', clear_path(path)+'OUTCAR')[0][-1].split()[5:8])
    except :
        return None


def get_Ef(PATH, same_dir = True,verbose = False, fileformat = 'OUTCAR'):
    tab = []
    if fileformat == 'OUTCAR' :
        try :
            if same_dir == False :
                for i in os.listdir(PATH+'..'):
                    if word in i :
                        tab.append(i)
                if len(tab) == 0 :
                    print('no accurate DOSCAR has been found, setting Ef from BS..., Ef =')
                    Ef = float(grep("fermi", PATH+"OUTCAR")[0][-1].split(":")[1].split()[0])
                    print(Ef)
                    return Ef
                elif len(tab) >= 2 :
                    print('many accurate DOSCAR has been found, setting Ef from :')
                    print(tab[0].split('/')[-1], ', Ef =')
                    Ef = float(grep("fermi", PATH+'../'+tab[-1]+'/'+"OUTCAR")[0][-1].split(":")[1].split()[0])
                    print(Ef)
                    return Ef
                else :
                    print('accurate DOSCAR has been found, setting Ef from :')
                    print(tab[0].split('/')[-1], ', Ef =')
                    Ef = float(grep("fermi", PATH+'../'+tab[-1]+'/'+"OUTCAR")[0][-1].split(":")[1].split()[0])
                    print(Ef)
                    return Ef
            else :
                #Ef = float(grep("fermi", PATH+"OUTCAR")[0][-1].split(":")[1].split()[0])
                Ef = float(grep("fermi", PATH+"OUTCAR")[0][-1].split(":")[1].split()[0])
                if verbose : print("Ef = %f" %(Ef))
                return Ef
        except :
            #row = get_row(6, PATH+'DOSCAR')
            #Ef = float(row.split()[-2])
            Ef = get_Ef_from_procar(PATH)
            return Ef
    elif fileformat == 'DOSCAR' :
        try :
            row = get_row(6, PATH+'DOSCAR')
            Ef = float(row.split()[-2])
            return Ef
        except :
            Ef = get_Ef_from_procar(PATH)
            return Ef
    else :
        Ef = get_Ef_from_procar(PATH)
        return Ef


def tail(f,n=1) :
    string = os.popen('tail -n '+str(n)+' '+'"'+f+'"').read()
    return string

def get_E0(path):
    pos = grep('sigma',path)[0][-1]
    E0 = float(pos.split()[-1])
    return E0

def get_direction(DIR):
    if DIR == 'x' : 
        return 0
    elif DIR == 'y':
        return 1
    elif DIR == 'z':
        return 2

def ERROR(message):
    try :
        raise ValueError(message)
    except ValueError:
        #print(message)
        sys.exit(message)

def get_IF(n_at,PLANE):
    c = 0
    for p in PLANE:
        if n_at in p :
            return c
        else :
            c+= 1
    return None

def get_ISPIN(path) :
    ISPIN = int(grep('ISPIN', clear_path(path)+'OUTCAR')[0][0].split()[2])
    return ISPIN

def get_spin(path, rtab = True):
    name = clear_path(path,space = False)+'OUTCAR'
    S = int(grep('ISPIN',name)[0][0].split()[2])
    try : 
        SOC = grep('LSORBIT',name)[0][0].split()[2]
        if SOC == 'F' :
            SOC = False
        else : 
            SOC = True
    except :
        SOC = False
    if rtab == True :
        if S == 1 :
            if SOC == True :
                return ['','_x','_y','_z']
            else :
                return ['']
        else :
            return ['_up','_dn']
    else :
        return S,SOC
        

