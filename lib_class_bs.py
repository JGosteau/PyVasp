#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 09:27:12 2019

@author: gosteau
"""

from lib_poscar import *

# Import C library to split PROCAR
path_program_split = "/home/gosteau/eclipse-workspace/PROCAR_split/"

HS_POINTS = {'G':[0.,0.,0.], 'X':[0.5,0.,0.], 'M':[0.5,0.5,0.], 'R':[0.5,0.0,0.5], 'X/3' : [0.165,0,0], 'M/3' : [0.165,0.165,0], 'A' : [0.5,0.5,0.5], 'Z' : [0,0,0.5]}

        
HS_ref = {'G' : [0.0,0.0,0.0], 
          'X' : [0.5,0.0,0.0], 'Y' : [0.0,0.5,0.0], 'Z' : [0.0,0.0,0.5],
          'M' : [0.5,0.5,0.0], 'N': [0.0,0.5,0.5], 'O': [0.5,0.0,0.5], 'R' : [0.5,0,0.5],
          'A' : [0.5,0.5,0.5]}
"""
HS_ref = {'G' : [0.0,0.0,0.0], 
          'A' : [0.0,0.0,0.5], 'I2' : [0.49949892,0.49949892,0.5], 'I' : [-0.49949892,0.50050108,0.5],
          'M2' : [-0.5,0.5,0.5], 'Y': [0.5,0.5,0], 'L2': [0.0,0.5,0.5], 'V2' : [0,0.5,0]}
"""

"""
Commande selection bande dans procar
def get_band():
    command="grep -A %d -E \"band +%d \" %s" %(n_atoms+1+space, band, PROCAR_name)

Commande obtenir projection:
def get_proj():
    command = "awk '$1 == \"tot\"{printf \"    0  \"; for(i=2; i<= NF; i++){printf $i\"  \"}print ""}$1 ~ /^[0-9]*$/ && NF != 0 {print $0}' "

"""


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
C = np.array([[255,125,0]])
cmap_o = mpl.colors.ListedColormap(C/255.0)
cmap_o.set_under('k',alpha=0)

def get_k(i) :
    if i%3 == 0:
        return 'kx'
    elif i%3 == 1:
        return 'kx'
    elif i%3 == 2:
        return 'kz'

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
    def __init__(self,name= None, kx=0.0, ky = 0.0,kz = 0.0, erase = True):
        try :
            sp = name.split('/')
            hs = HS(sp[0])/float(sp[1])
            self.__dict__ = hs.__dict__
        except :
            if name in HS_ref.keys() and erase:
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
"""
HS_points = {}
for i in HS_ref.keys() :
    hs = HS(i)
    HS_points[i] = hs
"""


def compile_split():
    ffi_procar_info = FFI()
    ffi_split = FFI()
    with open(path_program_c+"procar_info.h") as f :
        ffi_procar_info.cdef(f.read())
    ffi_procar_info.set_source("_procar_info_",
        """
        #include "procar_info.h"
        """, sources=["procar_info.c"])
    path_procar_info = ffi_procar_info.compile(tmpdir=path_program_c)
    lib_procar_info = ffi_procar_info.dlopen(path_procar_info)
    ffi_split.include(ffi_procar_info)
    with open(path_program_c+"split.h") as f :
        ffi_split.cdef(f.read())
    ffi_split.set_source("_split_",
        """
        #include "procar_info.h"
        #include "split.h"
        """, sources=["split.c"])
    path_split = ffi_split.compile(tmpdir=path_program_c)
    lib_split = ffi_split.dlopen(path_split)
    return lib_split, lib_procar_info

def get_lib_split():
    path_split = path_program_c + '_split_.so'
    path_procar_info = path_program_c + '_procar_info_.so'
    
    if os.path.isfile(path_split) and os.path.isfile(path_procar_info) :
        sys.path.append(path_program_c)
        from _split_ import lib as lib_split
        from _procar_info_ import lib as lib_procar_info
    else :
        lib_split, lib_procar_info = compile_split()
    return lib_split, lib_procar_info
#


def jobinfo(job):
    if type(job) == int :
        jobid = job
    else :
        jobid = int(job.split()[-1])
    sent = "scontrol show job %d" %(jobid)
    process = os.popen(sent)
    res = process.read()
    process.close()
    return res

def sbatch_bands(path, path_program_split, time = "00:10:00", N= 1, n = 1):
    string = """#!/bin/bash
#SBATCH -N %d
#SBATCH -n %d
#SBATCH --ntasks-per-node=%d
#SBATCH --ntasks-per-core=1
#SBATCH -J bands_split
#SBATCH --time=%s
#SBATCH --mail-user=none


#Positionnement de l'environnement
module purge
module load intel/18.2
module load intelmpi/18.2
#var environnements
export OMP_NUM_THREADS=1
export MALLOC_MMAP_MAX_=0
export MALLOC_TRIM_THRESHOLD_=-1
export FORT_BUFFERED=true
export decfort_dump_flag=true
export I_MPI_DAPL_SCALABLE_PROGRESS=1
export I_MPI_DAPL_TRANSLATION_CACHE=1
export MKL_CBWR=Auto

#repertoire d'installation des binaires VASP 5.4.4
export path='%s'
export SPLIT='%s'
cd $path
cp $SPLIT split.exe

#commande de lancement du binaire MPI (MPI pur, non multi-threade)
time srun $(placement 1 1 ) split.exe $path
echo 'Done'
rm split.exe
jobinfo $SLURM_JOBID
infoincidentjob
""" %(N,n*N,n,time, path, path_program_split+"split.exe")
    #print(string, string[573])
    with open(path+'bands.sh', 'w') as f :
        f.write(string)
    sent = "sbatch %sbands.sh" %(path)
    process = os.popen(sent)
    job = process.read()
    #print(job)
    process.close()
    return job



def create_bands(path, path_program_split, calmip = True,**kwargs):
    """
    Permet de découper le PROCAR par bandes dans le dossier BANDS, 
    il utilise le programme c "split.exe" situé dans le dossier path_program_split
    Exemple d'utilisation :
        path = "chemin vers le dossier où le PROCAR est present"
        create_bands(path)
    """
    path = clear_path(path)
    path_program_split = clear_path(path_program_split)
    if not os.path.isdir(path+'BANDS'):
        os.mkdir(path+'BANDS')
    else :
        shutil.rmtree(path+'BANDS')
        os.mkdir(path+'BANDS')
    if calmip == False :
        sent = path_program_split+"split.exe \""+path+"\""
        process = os.popen(sent)
        job = process.read()
        #print(job)
        process.close()
    else :
        job = sbatch_bands(path,path_program_split, **kwargs)
    return job
    
def cat_procar2(path, path_program_split, procar1, procar2, procarf = 'PROCAR'):
    path = clear_path(path)
    command = "grep -no k-points %s" %(path+procar1)
    process = os.popen(command)
    n1 = int(process.read().split()[-1].split(':')[0])
    process.close()
    command = "grep -no k-points %s" %(path+procar2)
    process = os.popen(command)
    n2 = int(process.read().split()[-1].split(':')[0])
    process.close()
    
    command = 'split -l%d %s %s' %(n1-2, path+procar1, path+'tmp_p1_')
    os.system(command)
    command = 'split -l%d %s %s' %(n2-2, path+procar2, path+'tmp_p2_')
    os.system(command)
    
    command = path_program_split+"cat_procar.sh %s%s %s%s %s%s" %(path, 'tmp_p1_aa', path, 'tmp_p2_aa', path, 'tmp_'+procarf+'_up')
    os.system(command)
    command = path_program_split+"cat_procar.sh %s%s %s%s %s%s" %(path, 'tmp_p1_ab', path, 'tmp_p2_ab', path, 'tmp_'+procarf+'_dn')
    os.system(command)
    
    command = "cat %s %s > %s" %(path+'tmp_'+procarf+'_up', path+'tmp_'+procarf+'_dn', path+procarf)
    os.system(command)
    command = "rm %stmp_p1_* %stmp_p2_* %stmp_%s*" %(path,path,path, procarf)
    os.system(command)    
    
def cat_procar(path, path_program_split, procar1, procar2, procarf = 'PROCAR', ISPIN = 1):
    if ISPIN == 1 :
        path = clear_path(path)
        command = path_program_split+"cat_procar.sh %s%s %s%s %s%s" %(path, procar1, path, procar2, path, procarf)
        os.system(command)
    else :
        cat_procar2(path, path_program_split, procar1, procar2, procarf)
   
    
    
    #command = path_program_split+"cat_procar.sh %s%s %s%s %s%s" %(path, procar1, path, procar2, path, procarf)
    #os.system(command)

path = "/home/gosteau/Documents/These/1ere année/STO/OUTPUT_dir_ncl/3.9/BS_SOC/"
n_ions = 5
n_kpts = 300
n_bands = 54

ion = 0
test_orbs  = np.array('s     py     pz     px    dxy    dyz    dz2    dxz  x2-y2    tot'.split())







def get_bands_limit(path, Elim, Ef = 0, spin = '',**kwargs):
    Emin, Emax = Elim
    #bands, occ = kwargs.get('bands', get_band_info(path, **kwargs))
    try : 
        bands, occ = get_band_info(path, spin = spin)
    except : 
        used_bands_1 = get_bands_limit(path,Elim,Ef, '_up')
        used_bands_2 = get_bands_limit(path,Elim,Ef, '_dn')
        used_bands = np.append(used_bands_1, used_bands_2)
        return np.arange(np.min(used_bands), np.max(used_bands)+1)
    used_bands = []
    i = 0
    while(np.min(bands[i])-Ef <= Emax):
        if(np.max(bands[i])-Ef >= Emin):
            used_bands.append(i+1)
        i += 1
    return np.array(used_bands)



class CONTRIB():      
    """
    La classe CONTRIB permet d'obtenir la contribution des atomes et des orbitales 
    d'une bande b en fonction des k-points.
    Exemple d'utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        band = b # n° de la bande b
        spin = '' # /'_up' / '_dn' / '_sx' / '_sy' / '_sz' # par défault la valeur est ''
        contrib = CONTRIB(path, band, spin)
        
    L'objet créé contient plusieurs informations :
        - contrib.band : correspond à la bande utilisée.
        - contrib.energy : Donne l'énergie de la bande en fonction des k-points.
        - contrib.occ : Donne l'occupation de la bande en fonction des k-points.
        - contrib.n_kpts : nombre de k-points
        - contrib.orbs : orbitales utilisées
        - contrib.n_ions : nombre d'ions
        - contrib.contrib : contribution des atomes et orbitales en fonction des k-points
    
    Pour obtenir la contribution de l'atome 'at'(1,2,3,...) et de l'orbitale 'orb' ('tot', 'px','py','pz',... ) :
        contrib[at][orb]
    Pour obtenir la contribution de l'atome 'at1' + 'at2' et de l'orbitale 'orb1' + 'orb2' ('tot', 'px','py','pz',... ) :
        contrib['at1+at2']['orb1+orb2'] # Attention 'at1+at2' doit être une chaine de caractère, idem pour 'orb1+orb2'
    Pour obtenir la contribution de l'atome 'at' et de l'orbitale 'orb' au point 'k' (0,1,2,...) :
        contrib[at][orb][k]
    Pour obtenir les contributions de tous les atomes et orbitales au point 'k' :
        contrib_k = contrib.print_tokpt(k)
    """
    
    def get_orb(self,orb):
        try :
            #return np.where(test_orbs == orb)[0][0]
            return np.where(self.orbs == orb)[0][0]
        except :
            string = orb +' not in [' 
            for l in self.orbs :
                string += l+', '
            string += ']'
            raise Exception(string)
    def test_list_orbs(self, liste):
        new = []
        if type(liste) == str :
            for l in liste.split('+') :
                new.append(self.get_orb(l))
        return new
    
    def test_list_ions(self, liste):
        new = []
        if type(liste) == str :
            if liste == 'tot' :
                new.append(0)
            else :
                for l in liste.split('+') :
                    new.append(int(l))
        elif type(liste) == int :
            new.append(liste)
        else :
            for i in liste :
                if i == 'tot' :
                    new.append(0)
                else :
                    new.append(int(i))
        return new
                
    
    def __init__(self, path, band, spin  ="", ions = None, info_procar = None):
        path = clear_path(path)
        if info_procar == None :
            info_procar = INFO_PROCAR(path)
        self.band = band
        def string_spin(spin):
            if 'x' in spin :
                return '_sx'
            elif 'y' in spin :
                return '_sy'
            elif 'z' in spin :
                return '_sz'
            elif 'up' in spin :
                return '_up'
            elif 'dn' in spin :
                return '_dn'
            else :
                return ''
        self.spin = string_spin(spin)
        tmp = get_band_info(path, spin = self.spin)
        self.energy = tmp[0][self.band-1]
        self.occ = tmp[1][self.band-1]
        self.n_kpts = info_procar.n_kpts
        #self.n_kpts = 280
        self.orbs = info_procar.orbs
        self.n_ions = info_procar.n_ions
        
        self.second_access = False
        
        if ions == None :
            ions = range(0,self.n_ions+1)
        ions = self.test_list_ions(ions)
        self.contrib = {}
        for i in ions :
            self.contrib[i] = []
        f = open(path+"BANDS/BAND_"+str(band)+spin, 'r')
        l = 1
        for line in f :
            for i in ions :
                if l%(self.n_ions+1) == i :
                    self.contrib[i].append(np.float64(line.split()[1:]))
            l += 1
        for i in ions :
            self.contrib[i] = np.transpose(np.array(self.contrib[i]))
        f.close()
    
    def tokpt(self, kpt, extra_orbs = True):
        """
        Permet d'obtenir toutes les contributions (atomes et orbitales) au point k :
            - kpt : point 'k' souhaité (entier)
            - extra_orbs : rajoute les contributions p et d. (boolean)
        """
        if extra_orbs:
            l = len(self.orbs)+2
        else :
            l = len(self.orbs)
        tab = np.zeros((self.n_ions+1, l))
        for i in range(self.n_ions+1):
            for j in range(len(self.orbs)):
                tab[i][j] = self.contrib[i][j][kpt]
                if extra_orbs :
                    if 'p' in self.orbs[j]:
                        tab[i][len(self.orbs)] += tab[i][j]
                    elif 'd' in self.orbs[j]:
                        tab[i][len(self.orbs)+1] += tab[i][j]
        return tab
    
    
    def print_tokpt(self, kpt,ions = None, extra_orbs = False, return_tab = True, return_string = False, kpoints = None, poscar = None, Ef = 0, supp_string = None):
        """
        Permet d'afficher toutes les contributions (atomes et orbitales) au point k :
            kpt : point 'k' souhaité (entier)
            extra_orbs : rajoute les contributions p et d. (boolean)
        """
        if Ef == None :
            Ef = 0
        string = ''
        tab = self.tokpt(kpt, extra_orbs = extra_orbs)
        if extra_orbs == True :
            sup = 2
        else :
            sup = 0
        lmax = int(4 + 7*(len(self.orbs)+sup))
        string += '-'*lmax+'\n'
        
        #l_num = 15+11
        sentence = ' K-POINT %4d - BAND %4d%s ' %(kpt, self.band,self.spin)
        l_num = len(sentence)
        l = int((lmax-l_num)/2)
        string += '-'*l
        string += sentence
        string += '-'*(lmax-l-l_num) +"\n"
        
        if kpoints != None :
            #l = int((4 + 7*(len(self.orbs)+sup)-3-3*10 )/2)
            sentence = ' %1.8f %1.8f %1.8f ' %(kpoints['kx'][kpt],kpoints['ky'][kpt],kpoints['kz'][kpt])
            l_num = len(sentence)
            l = int((lmax-l_num)/2)
            string += '-'*l
            string += sentence
            string += '-'*(lmax-l-l_num) +"\n"
        if supp_string != None :
            l_num = len(supp_string)
            l = int((lmax-l_num)/2)
            string += '-'*l
            string += supp_string
            string += '-'*(lmax-l-l_num) +"\n"
        #l = int((4 + 7*(len(self.orbs)+sup)-25-15 )/2)
        sentence = ' Energy %3.8f eV - Occupancy %3.3f ' %(self.energy[kpt]-Ef,self.occ[kpt])
        l_num = len(sentence)
        l = int((lmax-l_num)/2)
        string += '-'*l
        string += sentence
        string += '-'*(lmax-l-l_num) +"\n"
        
        string += '-'*lmax +"\n"
        
        string += " ion"
        #print("ion",end = '')
        for orb in self.orbs :
            string += "  %5s" %(orb)
            #print("  %5s" %(orb),end = '')
        if extra_orbs:
            string += "  %5s  %5s" %('p', 'd')
            #print("  %5s  %5s" %('p', 'd'),end = '')
        string += "\n"
        #print("")
        for i in range(1,self.n_ions+1):
            if type(poscar) == type(None) :
                string += "%4d" %i
            else :
                string += "%4s" %(poscar.get_count_at(i))
            #print("%3d" %i, end='')
            for j in range(len(self.orbs)):
                if tab[i][j] < 0 :
                    string += ' -%1.3f' %(abs(tab[i][j]))
                else : 
                    string += "  %1.3f" %(abs(tab[i][j]))
                #print("  %1.3f" %(tab[i][j]),end = '')
            if extra_orbs:
                string += "  %1.3f  %1.3f" %(tab[i][len(self.orbs)],tab[i][len(self.orbs)+1])
                #print("  %1.3f  %1.3f" %(tab[i][len(self.orbs)],tab[i][len(self.orbs)+1]),end = '')
            string += "\n"
            #print("")
        string += " tot"
        #print("tot", end = '')
        for j in range(len(self.orbs)):
            
            if extra_orbs:
                if 'p' in self.orbs[j]:
                    tab[i][len(self.orbs)] += tab[i][j]
                elif 'd' in self.orbs[j]:
                    tab[i][len(self.orbs)+1] += tab[i][j]
            if tab[0][j] < 0 :
                string += ' -%1.3f' %(abs(tab[0][j]))
            else : 
                string += "  %1.3f" %(abs(tab[0][j]))
            
            #string += "  %1.3f" %(tab[0][j])
            #print("  %1.3f" %(tab[0][j]),end = '')
        if extra_orbs:
            string += "  %1.3f  %1.3f" %(tab[0][len(self.orbs)],tab[0][len(self.orbs)+1])
            #print("  %1.3f  %1.3f" %(tab[0][len(self.orbs)],tab[0][len(self.orbs)+1]),end = '')
        string += "\n"
        #print("")
        if return_string and return_tab :
            return tab,string
        elif return_string :
            return string
        elif return_tab :
            print(string)
            return tab
        else :
            print(string)
    
    def __getitem__(self, key):
        if not self.second_access:
            self.first_key = key
            self.second_access = True
            return self
        else:
            self.second_access = False
            ions = self.first_key
            orbs = key
            ions = self.test_list_ions(ions)
            orbs = self.test_list_orbs(orbs)
            if len(ions) == 1 and len(orbs) == 1 :
                return self.contrib[ions[0]][orbs[0]]
            else :
                result = np.zeros(self.n_kpts)
                for i in ions:
                    for l in orbs :
                        result += self.contrib[i][l]
            return result
      
    def __repr__(self) :
        string = ''
        for kpt in range(self.n_kpts):
            string += self.print_tokpt(kpt, return_string = True, return_tab = False)
        return string        

        
class KPOINT():
    def __init__(self, kpt=None, kx=None,ky=None,kz=None,weight=None):
        self.kpt = kpt
        self.kx = kx
        self.ky = ky
        self.kz = kz
        try :
            self.norm = (kx**2+ky**2+kz**2)**(1/2.)
        except :
            self.norm = None
        self.weight = weight
        
def new_ticks(kpoints, index, knorm, HSs = HS_points):
    ticks, label = kpoints.search_all_hs(index, HSs = HSs)
    new_tick = []
    new_label = []
    uniques = np.unique(knorm[ticks], return_counts=True)
    for i in np.arange(len(uniques[1])):
        c = uniques[1][i]
        u = uniques[0][i]
        if c == 1 :
            j = np.where(knorm[ticks] == knorm[ticks][i])[0][0]
            new_tick.append(ticks[j])
            new_label.append(label[j])
        if c == 2 :
            wh = np.where(knorm[ticks] == u)[0]
            #print(i,u, wh, ticks[wh[0]], knorm[ticks][i])
            if ticks[wh[0]] <  ticks[wh[1]] :
                text = label[wh[0]]+'/'+label[wh[1]]
            else :
                text = label[wh[1]]+'/'+label[wh[0]]
            #print(text)
            new_tick.append(ticks[wh[0]])
            new_label.append(text)
    return new_tick, new_label

class KPOINTS():
    """
    Classe permettant d'obtenir les coordonnées des différents k-points.
    Exemple d'utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        ff = 'PROCAR' # default = 'PROCAR' / 'OUTCAR' / 'wannier'
        kpoints = KPOINTS(path, fileformat = ff)
        
        kpoints # Afficher tous les k-points
        
        kpoints.loc[kpoints['kz'] == 0] # Affiche tous les k-points qui respectent kz = 0
        kpoints.loc[kpoints['kx'] == kpoints['ky']] # Affiche tous les k-points qui respectent kx = ky
        kpoints.loc[kpoints['kx'] == kpoints['ky']].index # Donne l'index des k-points qui respectent kx = ky
        
        kpoints.loc[range(0,10)] # Affiche les 10 premiers k-points
        
        kpoints['kx'] # Affiche la coordonnée kx de tous les k-points
        kpoints['kx'][kpoints['kz'] == 0] # Affiche la coordonnée kx des k-points respectant kz = 0
        kpoints['kx'].loc[kpoints['kz'] == 0] # Affiche la coordonnée kx des k-points respectant kz = 0
        kpoints.loc[kpoints['kz'] == 0]['kx'] # Affiche la coordonnée kx des k-points respectant kz = 0
        kpoints[['kx','ky']] # Affiche les coordonnées kx et ky de tous les k-points
        
        kpoints.search_hs(HS_pt) # cherche la position du point de haute symétrie HS_pt
        kpoints.remove_double() # Donne l'index des kpoints en supprimant les k-points qui se suivent et qui ont les mêmes coordonnées.
        
        kpoints.search_all_hs() # renvoie tous les points de HS
    """
    
    def load_outcar(self,path):
        if self.reciprocal :
            #l_kpt=int(grep("k-points in reciprocal lattice and weights: K", path + "/OUTCAR")[1][0])
            l_kpt=int(grep("k-points in reciprocal lattice and weight", path + "/OUTCAR")[1][0])
            #l_kpt=int(grep("Following reciprocal coordinates:", path + "/OUTCAR")[1][0])+1
        else : 
            l_kpt=int(grep("k-points in units of 2pi/SCALE and weight:", path + "/OUTCAR")[1][0])
        kpt = pd.read_csv(path + "OUTCAR", skiprows=l_kpt, nrows=self.nkpt, names=['kx','ky','kz','w_k'], delim_whitespace=True)
        kpt['kpt'] = pd.Series(range(1,self.nkpt+1),index=kpt.index)
        return kpt
    
    def __init__(self, path=None, reciprocal = True, fileformat = 'PROCAR'):
        if path != None :
            if fileformat == 'OUTCAR' :
                self.nkpt = int(grep("NKPTS" , path + "OUTCAR")[0][0].split()[3]) 
                self.reciprocal = reciprocal
                self.data = self.load_outcar(path)
            elif fileformat == 'PROCAR' :
                self.data = pd.read_csv(path+"BANDS/kpts_tmp", names = ["kx","ky","kz","weight"], delim_whitespace = True)
                self.nkpt = len(self.data)
            elif fileformat == 'wannier' :
                try :
                    self.nkpt = int(get_row(1, path+'wannier90_band.kpt').split()[0])
                    self.data = pd.read_csv(path+"wannier90_band.kpt", names = ["kx","ky","kz","weight"], delim_whitespace = True, skiprows = 1, nrows = self.nkpt)
                except :
                    self.nkpt = int(get_row(1, path+'wannier90.up_band.kpt').split()[0])
                    self.data = pd.read_csv(path+"wannier90.up_band.kpt", names = ["kx","ky","kz","weight"], delim_whitespace = True, skiprows = 1, nrows = self.nkpt)
                
            self.data['k'] = np.sqrt(self['kx']**2+self['ky']**2+self['kz']**2)
            
        else :
            self.nkpt = None
            self.reciprocal = reciprocal
            self.data = pd.DataFrame(columns=['kpt','kx','ky','kz','w_k'])    
        self.loc = self.data.loc
        self.u, self.v , self.w = self.get_uvw(path)

    def get_uvw(self,path):
        try :
            file = clear_path(path) + 'OUTCAR'
            n = grep("reciprocal lattice vectors", file)[1][0]
            u = np.float64(np.array(get_row(n+1, file).split())[3:])*np.pi*2
            v = np.float64(np.array(get_row(n+2, file).split())[3:])*np.pi*2
            w = np.float64(np.array(get_row(n+3, file).split())[3:])*np.pi*2
        except :
            u = np.array([1,0,0])
            v = np.array([0,1,0])
            w = np.array([0,0,1])
        return np.array([u,v,w])
 
    def get_knorm(self, kpt_line,u = None,v = None,w = None):
        if type(u) == type(None) :
            u = self.u
        if type(v) == type(None) :
            v = self.v
        if type(w) == type(None) :
            w = self.w
        
        KNORM = np.array([])
        KPTS = []
        i = 0
        for line in kpt_line:
            knorm_line = []
            kpt = line[0]
            KPTS.append(kpt)
            kx = self.loc[kpt]['kx']*u
            ky = self.loc[kpt]['ky']*v
            kz = self.loc[kpt]['kz']*w
            k_ref = kx+ky+kz
            knorm_line.append(0)
            
            ini = line[0]
            end = line[1]
            if line[0] > line[1] :
                p= 0
                f = 0
                incr = -1
            else :
                incr = 1
                p = 1
                f = 1
            
            #KPTS += range(ini+p,end+f, incr)
            for kpt in range(ini+p,end+f, incr):
                KPTS.append(kpt)
                #print(line)
                kx = self.loc[kpt]['kx']*u
                ky = self.loc[kpt]['ky']*v
                kz = self.loc[kpt]['kz']*w
                k = kx+ky+kz
                knorm_line.append(np.linalg.norm(k-k_ref))
            if i == 0 :
                KNORM = np.append(KNORM,knorm_line)
                i = 1
            else :
                KNORM = np.append(KNORM,knorm_line+KNORM[-1])
        return KPTS, KNORM
    
    def search_hs(self, hs, unique = False):
        """
        Retourne l'index du point de Haute Symétrie hs dans la grille des k-points.
        
        Utilisation :
            hs = HS_pt('G')
            kpoints.search_hs(hs)
        """
        points = []
        if type(hs) == str or type(hs) == dict :
            hs = HS_pt(hs)
            loc = np.array(self.data.loc[self.data['kx'] == hs.kx]
                .loc[self.data['ky'] == hs.ky]
                .loc[self.data['kz'] == hs.kz].index)
            if unique :
                if len(loc) >0:
                    return loc[0]
                else :
                    return loc
            else :
                return np.array(self.data.loc[self.data['kx'] == hs.kx]
                .loc[self.data['ky'] == hs.ky]
                .loc[self.data['kz'] == hs.kz].index)
        elif type(hs) == list or type(hs) == np.ndarray :
            for h in hs :
                points.append(self.search_hs(h))
            return points
        elif type(hs) == type(None) :
            return []
        else :
            loc = np.array(self.data.loc[self.data['kx'] == hs.kx]
                .loc[self.data['ky'] == hs.ky]
                .loc[self.data['kz'] == hs.kz].index)
            if unique :
                if len(loc) >0:
                    return loc[0]
                else :
                    return None
            else :
                return loc
      
    def search_all_hs(self, index = None):
        ticks = np.array([])
        label = np.array([])
        if type(index) == type(None) :
            index = self.data.index
        for i in HS_points.keys() :
            hs = HS_points[i]
            hs_pos_in_kpt = self.search_hs(hs)
            k = 0
            while k < len(hs_pos_in_kpt) :
                kpt = hs_pos_in_kpt[k]
                if kpt not in index :
                    hs_pos_in_kpt = np.delete(hs_pos_in_kpt,k)
                else :
                    k += 1
            #ax.vlines(hs_pos_in_kpt, *ylim, **kwargs)
            test = []
            #print(hs_pos_in_kpt)
            for j in hs_pos_in_kpt :
                #print(j, np.where(kpt == j))
                test.append(np.where(index == j)[0][0])
            test = np.array(test)
            #ticks = np.append(ticks, hs_pos_in_kpt)
            ticks = np.append(ticks, test)
            label = np.append(label, [i]*len(hs_pos_in_kpt))
        return ticks, label       

    #def search_ticks(index, knorm):
    def search_all_hs(self, index = None, knorm = None, HSs = HS_points):
        if type(index) == type(None) :
            index = self.data.index
        ticks = np.array([])
        label = np.array([])
        for i in HSs.keys() :
            hs = HSs[i]
            hs_pos_in_kpt = self.search_hs(hs)
            k = 0
            while k < len(hs_pos_in_kpt) :
                kpt = hs_pos_in_kpt[k]
                if kpt not in index :
                    hs_pos_in_kpt = np.delete(hs_pos_in_kpt,k)
                else :
                    k += 1
            for j in hs_pos_in_kpt :
                wh = np.where(index == j)[0]
                ticks = np.append(ticks, wh)
                label = np.append(label, [i] * len(wh))
        return np.int0(ticks), label
    
    def ticks_for_plot(self, index, knorm):
        """
        Renvoi les points de HS dans un format adapté au plot
        """
        ticks, label = self.search_all_hs(index)
        new_tick = []
        new_label = []
        uniques = np.unique(knorm[ticks], return_counts=True)
        for i in np.arange(len(uniques[1])):
            c = uniques[1][i]
            u = uniques[0][i]
            if c == 1 :
                j = np.where(knorm[ticks] == knorm[ticks][i])[0][0]
                new_tick.append(ticks[j])
                new_label.append(label[j])
            if c == 2 :
                wh = np.where(knorm[ticks] == u)[0]
                #print(i,u, wh, ticks[wh[0]], knorm[ticks][i])
                if ticks[wh[0]] <  ticks[wh[1]] :
                    text = label[wh[0]]+'/'+label[wh[1]]
                else :
                    text = label[wh[1]]+'/'+label[wh[0]]
                #print(text)
                new_tick.append(ticks[wh[0]])
                new_label.append(text)
        return new_tick, new_label    
    
    def auto_lines(self):
        """
        Determine automatiquement les points de hautes symétries
        """
        index, hs = self.search_all_hs()
        hs = hs[np.argsort(index)]
        hs_lines = [] 
        for i in range(0,len(hs),2):
            hs_lines.append([hs[i],hs[i+1]])
        return hs_lines
    
    def __getitem__(self, index):
        return self.data[index]
    
    def __len__(self):
        return self.nkpt
    
    def __repr__(self):
        return "%s" %str(self.data)
           
    def remove_double(self, index, knorm):
        """
        Renvoie l'index des k-points en supprimant les k-points qui se suivent et qui ont 
        les mêmes coordonnées.
        """
        i = 0
        li = len(index)
        while i < li-1 :
            kx1, kx2 = self['kx'][index[[i,i+1]]]
            ky1, ky2 = self['ky'][index[[i,i+1]]]
            kz1, kz2 = self['kz'][index[[i,i+1]]]
            if kx1 == kx2 and ky1 == ky2 and kz1 == kz2 :
                knorm = np.delete(knorm, i+1)
                index = np.delete(index, i+1)
                li -= 1
            else :
                i += 1
        return np.int0(index), knorm
    
    def get_one_dir(self,hs1, hs2, prev_knorm = 0, fax= 1 , fay = 1, faz = 1):
        """
        Renvoie l'index des points k ainsi que la norme 
        satisfaisant la direction de HS déterminé par hs1 et hs2
        """
        def sign(x):
            if x >= 0 :
                return True
            else :
                return False
            
        if type(hs1) != HS :
            hs1 = HS(hs1)
        if type(hs2) != HS :
            hs2 = HS(hs2)
        if hs2 == hs1 :
            return np.array([]),np.array([])
        hs = hs2-hs1
        count = 0
        not_zeros_dir = []
        cond = []
        for i in ['kx','ky','kz']:
            if hs[i] == 0 :
                cond.append(self[i] == hs2[i])
                count += 1
            else :
                not_zeros_dir.append(i)
        if count == 0 :
            cond.append(hs['kx'] * self['kx'] == hs['ky'] * self['ky'])
            cond.append(hs['kx'] * self['kx'] == hs['kz'] * self['kz'])
                
        elif count == 1 :
            dir1 = not_zeros_dir[0]
            dir2 = not_zeros_dir[1]
            cond.append(hs[dir1] * self[dir1] == hs[dir2] * self[dir2])
        #print count
        kpt = np.array(self.loc[cond[0]].loc[cond[1]].sort_values(['kx','ky','kz'], 
                       ascending = [sign(hs[i]) for i in ['kx','ky','kz']]).index)
        
        
        A = {'kx' : fax, 'ky': fay, 'kz' : faz}
        knorm = np.zeros(len(kpt)) 
        for i in ['kx','ky','kz'] : 
            knorm += ((self[i][kpt] - hs1[i])/A[i])**2
        knorm = np.sqrt(knorm) + prev_knorm
        return kpt, knorm
    
    def get_dirs(self, HLines, fax = 1, fay = 1, faz = 1):
        index = np.array([])
        knorm = np.array([])
        prev_knorm = 0
        for i in range(len(HLines)-1) :
            hs1 = HLines[i]
            if type(hs1) == str :
                hs1 = HS(hs1)
            hs2 = HLines[i+1]
            if type(hs2) == str :
                hs2 = HS(hs2)
            new_index, new_knorm = self.get_one_dir(hs1,hs2, prev_knorm = prev_knorm, fax=fax,fay=fay,faz=faz)
            index = np.append(index, new_index)
            knorm = np.append(knorm, new_knorm)
            prev_knorm = knorm[-1]
        return self.remove_double(np.int0(index), knorm)
        #return np.int0(index), knorm
        
    def get_dirs(self, HLines, fax = 1, fay = 1, faz = 1):
        index = np.array([])
        knorm = np.array([])
        prev_knorm = 0
        for i in range(len(HLines)) :
            hs1 = HLines[i][0]
            hs2 = HLines[i][1]
            
            new_index, new_knorm = self.get_one_dir(hs1,hs2, prev_knorm = prev_knorm, fax=fax,fay=fay,faz=faz)
            index = np.append(index, new_index)
            knorm = np.append(knorm, new_knorm)
            prev_knorm = knorm[-1]
        return self.remove_double(np.int0(index), knorm)
        #return np.int0(index), knorm
    
            

def onpick_base(event):
    thisline = event.artist
    xdata = thisline.get_xdata()
    ydata = thisline.get_ydata()
    ldata = thisline.get_label()
    ind = event.ind
    for i in range(len(xdata[ind])):
        x = xdata[ind][i]
        y = ydata[ind][i]
        print('label: ', ldata, ' ; x: ', x, ' ; y: ', y, ' ; ind :', i)
            
        



def rashba_calc_fit(path, b1= None, b2= None,R= None,  
                    ref_kpt = None, unit = 'normal',
                    verbose = True, 
                    PLOT = True,  figname = '' ,
                    bands_range = None):
    """
    Renvoie la valeur des fits linéaire, cubique et linéaire+cubique de la différence d'énergie
    entre deux bandes b1 et b2 dans une gamme de k-points R.
    - b1 : Première bande
    - b1 : Deuxième bande
    - R : index des k-points à utiliser pour les fits.
    
    - ref_kpt : valeur de référence pour le fit ie k-ref_kpt, format :
    ref_kpt = {'kx' : _, 'ky' : _, 'kz' : _}
    Par défault ref_kpt est donné par le premier point de R.
    - unit : unité pour le fit, 'normal' pour eV.A, 'reciprocal' pour eV/(2pi/a)
    
    - verbose : Donne plus de détails sur les fit si True
    
    - PLOT : affiche la dispersion des bandes ainsi que les fits de la différence d'énergie si True
    - figname : nom de la figure
    
    - bands_range = {'bands' : [[b1,b2], [b3,b4], ...], 'R' : [R1,R2]} : permet de faire les fits 
    sur plusieurs bandes si des bandes se croisent, (pas encore terminé).
    """
    kpoints = KPOINTS(path, fileformat = 'PROCAR')
    if unit != 'reciprocal' :
        poscar = POSCAR(path)
        ax, ay, az = poscar.diag_lattice
    else :
        ax = 1
        ay = 1
        az = 1
    def get_bands_and_kpt(b1,b2,R, ref_kpt = None):
        band_1 = get_band_info(path)[0][b1-1][R]
        band_2 = get_band_info(path)[0][b2-1][R]
        if type(ref_kpt) == type(None):
            ref_kpt = kpoints.loc[R[0]]
        new_kpt = kpoints.loc[R]
        new_kpt['k'] = np.sqrt(((new_kpt['kx']-ref_kpt['kx'])*2*np.pi/ax)**2 + ((new_kpt['ky']-ref_kpt['ky'])*2*np.pi/ay)**2 + ((new_kpt['kz']-ref_kpt['kz'])*2*np.pi/az)**2)
        return band_1, band_2, new_kpt
    
    
    if b1 != None :
        band_1, band_2, new_kpt = get_bands_and_kpt(b1,b2,R, ref_kpt)
        #print(new_kpt)
    else :
        band_1 = np.array([])
        band_2 = np.array([])
        list_band = bands_range['bands']
        list_range = bands_range['R']
        ref_kpt = kpoints.loc[list_range[0][0]]
        new_kpt = pd.DataFrame(columns = ref_kpt.keys())
        for i in range(len(list_band)) :
            b1 = list_band[i][0]
            b2 = list_band[i][1]
            R = list_range[i]
            tmp_1, tmp_2, tmp_kpt = get_bands_and_kpt(b1,b2,R,ref_kpt)
            band_1 = np.append(band_1, tmp_1)
            band_2 = np.append(band_2, tmp_2)
            new_kpt = new_kpt.append(tmp_kpt)
        
        
    dE = band_2 - band_1
    
    # ---- FIT ----
        # 1) Function to fit :
    lin_fit = lambda x,a,c : a*x+c
    cub_fit = lambda x,a,c : a*x**3+c
    car_fit = lambda x,a,c : a*x**2+c
    mix_fit = lambda x,a,b,d,c : a*x**3+b*x**2+d*x+c
    mix2_fit = lambda x,a,b,c : a*x**3+b*x+c
    
        # 2) Apply the fit :
    popt_lin, _ = curve_fit(lin_fit, new_kpt['k'],dE)
    perr_lin = np.sqrt(np.diag(_))*100
    popt_cub, _ = curve_fit(cub_fit, new_kpt['k'],dE)
    perr_cub = np.sqrt(np.diag(_))*100
    popt_car, _ = curve_fit(car_fit, new_kpt['k'],dE)
    perr_car = np.sqrt(np.diag(_))*100
    popt_mix, _ = curve_fit(mix_fit, new_kpt['k'],dE)
    perr_mix = np.sqrt(np.diag(_))*100
    popt_mix2, _ = curve_fit(mix2_fit, new_kpt['k'],dE)
    perr_mix2 = np.sqrt(np.diag(_))*100
    
        # 3) Print result :
    if verbose :
        string = 'lin : a*k + b : \n  -a = %.5f eV/A (%.2f %%)\n  -b = %.5f eV (%.2f %%)\n' %(popt_lin[0],perr_lin[0],popt_lin[1],perr_lin[1])
        print(string)
        
        string = 'cub : a*k**3 + b : \n  -a = %.5f eV/A**3 (%.2f %%)\n  -b = %.5f eV (%.2f %%)\n' %(popt_cub[0],perr_cub[0],popt_cub[1],perr_car[1])
        print(string)
        
        string = 'mix : a*k**3 + b*k + c : \n  -a = %.5f eV/A**3 (%.2f %%)\n  -b = %.5f eV/A (%.2f %%)\n  -c = %.5f eV (%.2f %%)\n' %(popt_mix2[0],perr_mix2[0],popt_mix2[1],perr_mix2[1],popt_mix2[2],perr_mix2[2])
        print(string)
    
    # --- DISPLAY ---
    if PLOT == True:
        picker_res = 2
        fig = plt.figure(figname)
        ax = plt.subplot(1,2,1)
        ax.plot(new_kpt['k'], band_1, label = str(b1), color = 'blue', marker = '+', picker = picker_res)
        ax.plot(new_kpt['k'], band_2, label = str(b2), color = 'red', marker = '+', picker = picker_res)
        plt.legend()
        
        ax = plt.subplot(1,2,2)
        ax.plot(new_kpt['k'], dE, label = 'dE', color = 'black', marker = '+', picker = picker_res)
        plt.plot(new_kpt['k'], lin_fit(new_kpt['k'],*popt_lin), label = 'lin : a=%.3f'%(popt_lin[0]), linestyle = ':', color = 'orange', picker = picker_res)
        plt.plot(new_kpt['k'], cub_fit(new_kpt['k'],*popt_cub), label = 'cub : a=%.3f'%(popt_cub[0]), linestyle = ':', color = 'purple', picker = picker_res)
        plt.plot(new_kpt['k'], mix2_fit(new_kpt['k'],*popt_mix2), label = 'cub+lin : a=%.3f , b = %.3f'%(popt_mix2[0], popt_mix2[1]), linestyle = ':', color = 'black', picker = picker_res)
        plt.legend()
        
        fig.canvas.mpl_connect('pick_event', onpick_base) 
    return [popt_lin, popt_cub, popt_mix2], [perr_lin, perr_cub, perr_mix2]
     
                  
def Band_Structure_interp(path, Elim, Ef= None, kpoints = None, kpt = None, HS_dir = None, interp_kind = None,n_interp = 1000, orb = 'tot', atom = 'tot', spin = '',cmap = plt.get_cmap('seismic'), norm = mpl.colors.Normalize(0,1), k_norm = False, BANDS = None,picker = 1, Disp_BANDS = None, knorm = None, k_unit = 'rec'):
    """
    (Vielle version, il vaut mieux utiliser celle défini dans lib_display_BS)
    
    
    Permet de tracer la stucture de bande et faire des projections sur les atomes et orbitales.
    Actuellement le programme trace en fonction de l'index des kpoints et non la norme
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi (si None Ef est déterminé par le OUTCAR ou le DOSCAR du path)
    - kpoints : grille de kpoints si déjà calculé 
    - kpt : index des k-points à afficher
    - HS_dir : pas encore implémenté
    - interp_kind : type d'interpolation si None --> pas d'interpolation ('cubic' / 'linear' / None)
    - n_interp : nombre de point pour l'interpolation.
    - orb = 'tot' : orbitale projetée ('tot', 'px', 'px+py', ...)
    - at = 'tot' : atome projeté ('tot', 1, '1', '1+2', ...)
    - spin = '' : projection sur les spins '' / '_sx' / '_sy' / '_sz'
    - cmap : color_map utilisée
    - norm : norm pour la projection.
    
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Band_Structure_interp(path, Elim = [Emin, Emax])
    """
    def string_spin(spin):
        if 'x' in spin :
            return '_sx'
        elif 'y' in spin :
            return '_sy'
        elif 'z' in spin :
            return '_sz'
        else :
            return ''
    
    if Ef == None :
        Ef = get_Ef(path,same_dir=True)
    if type(Disp_BANDS) == type(None) :
        Disp_BANDS = get_bands_limit(path, Elim, Ef)
    if type(kpoints) == type(None):
        kpoints = KPOINTS(path, from_procar = True)
    if HS_dir == None :
        if type(kpt) == type(None):
            kpt = kpoints.remove_double(kpoints.data.index)
    #else :
    #    kpt = kpoints.get_dirs(HS_dir)
    edgecolor = 'none'
    size = 1
    marker = '.'
    if BANDS == None :
        BANDS = {}
    ref_glob = time.time()
    for b in Disp_BANDS:
        if b not in BANDS.keys() :
            ref = time.time()
            bands = CONTRIB(path, b, string_spin(spin))
            BANDS[b] = bands
        else : 
            bands = BANDS[b]
        
        #k = kpoints.loc[kpt].reset_index().index
        k = range(len(kpt))
        #k = kpt
        if interp_kind != None :
            new_bands = inter.interp1d(kpt,bands.energy[kpt], kind = interp_kind)
            new_orb = inter.interp1d(kpt,bands[atom][orb][kpt], kind = interp_kind)
            k_new = np.linspace(kpt[0],kpt[-1],n_interp)
            band_out = new_bands(k_new)
            orb_out = new_orb(k_new)
            k = k_new
        else :
            band_out = bands.energy[kpt]
            
            orb_out = bands[atom][orb][kpt]
        
        ref = time.time()
        points = np.array([k,band_out-Ef]).T.reshape(-1,1,2)
        segments = np.concatenate([points[:-1], points[1:]], axis = 1)
        lc = LineCollection(segments,cmap=cmap,norm = norm)
        lc.set_array(orb_out)
        plt.gca().add_collection(lc)
        plt.plot(k,band_out-Ef,c = 'black', marker = marker, alpha = 0, linewidth = 1, label = str(b), picker = picker)
        plt.xlim(k[0],k[-1])
        plt.ylim(*Elim)
    return BANDS


def search_sym_lines_and_plot(ax, kpoints, index,knorm,k_unit = 'rec', **kwargs):
    ylim = ax.get_ylim()
    ticks = np.array([])
    label = np.array([])
    for i in HS_points.keys() :
        hs = HS_points[i]
        hs_pos_in_kpt = kpoints.search_hs(hs)
        k = 0
        while k < len(hs_pos_in_kpt) :
            kpt = hs_pos_in_kpt[k]
            if kpt not in index :
                hs_pos_in_kpt = np.delete(hs_pos_in_kpt,k)
            else :
                k += 1
        #ax.vlines(hs_pos_in_kpt, *ylim, **kwargs)
        test = []
        #print(hs_pos_in_kpt)
        for j in hs_pos_in_kpt :
            #print(j, np.where(kpt == j))
            test.append(np.where(index == j)[0][0])
        test = np.array(test)
        #ticks = np.append(ticks, hs_pos_in_kpt)
        ticks = np.append(ticks, test)
        label = np.append(label, [i]*len(hs_pos_in_kpt))
    ax.set_xticks(ticks)
    ax.vlines(ticks, *ylim, **kwargs)
    ax.set_xticklabels(label)
    return ticks, label

    
            
def axes_BS_interp(path, Elim , Ef = None,at = 'tot', orb = 'tot', sp = '', cmap = plt.get_cmap('seismic'), path_Ef = None, HS_dir = None, interp_kind = None , ax = None, norm = None, new_fig = True, figsize = None, figname = "BS", BANDS = None, plot_ticks = True, picker = 1, kpoints = None, kpt = None, Disp_BANDS = None, poscar = None, fig = None, k_unit = 'rec'):
    """
    (Vielle version, il vaut mieux utiliser celle défini dans lib_display_BS : plot_BS_contrib)
    
    Permet de tracer la stucture de bande et faire des projections sur les atomes et orbitales
    en rajoutant les points de HS et le niveau de Fermi.
    Actuellement le programme trace en fonction de l'index des kpoints et non la norme
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi (si None Ef est déterminé par le OUTCAR ou le DOSCAR du path)
    - HS_dir : pas encore implémenté
    - interp_kind : type d'interpolation si None --> pas d'interpolation ('cubic' / 'linear' / None)
    - n_interp : nombre de point pour l'interpolation.
    - orb = 'tot' : orbitale projetée ('tot', 'px', 'px+py', ...)
    - at = 'tot' : atome projeté ('tot', 1, '1', '1+2', ...)
    - sp = '' : projection sur les spins '' / '_sx' / '_sy' / '_sz'
    - cmap : color_map utilisée
    - norm : norm pour la projection.
    - kpoints : grille de kpoints si déjà calculé 
    - kpt : index des k-points à afficher
    - figname : Nom de la figure.
    
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        axes_BS_interp(path, Elim = [Emin, Emax])
    """
    if kpoints == None :
        kpoints = KPOINTS(path, from_procar = True)
    if HS_dir == None :
        if type(kpt) == type(None) :
            kpt = kpoints.remove_double(kpoints.data.index)
    knorm = kpt
    if fig == None :
        if new_fig :
            fig = plt.figure(figname, figsize = figsize)
        
    def onpick(event):
        thisline = event.artist
        #thisline = event
        
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ldata = thisline.get_label()
        cdata = thisline.get_color()
        ind = event.ind
        k_ind = np.int0(xdata[ind])
        energy = np.float64(ydata[ind])
        band = np.int0(ldata)
        #print(band)
        #print("######## POINTS ########") 
        for i in range(len(k_ind)) :
            index = kpt[k_ind[i]]
            kx = kpoints['kx'][index]
            ky = float(kpoints['ky'][index])
            kz = float(kpoints['kz'][index])
            k = float(kpoints['k'][index])
            if len(k_ind) > 1:
                #e = energy[i]
                #print('point k(%d) = [%f ; %f ; %f ; %f]\nBand : %d\nEnergy = %f eV\n-----------' %(index,kx,ky,kz,k,band,e))
                BANDS[band].print_tokpt(index, kpoints = kpoints,poscar = poscar, Ef = Ef)
            else :
                #print('point k(%d) = [%f ; %f ; %f ; %f]\nBand : %d\nEnergy = %f eV' %(index,kx,ky,kz,k,band,energy))
                BANDS[band].print_tokpt(index, kpoints = kpoints, poscar = poscar, Ef = Ef)
    if new_fig :
        fig.canvas.mpl_connect('pick_event', onpick) 
        
    if path_Ef == None:
        path_Ef = path
    if Ef == None :
        Ef_SCF = get_Ef(path_Ef, same_dir = True)
    else :
        Ef_SCF = Ef
    
    #if HS_dir != None:
    #    kpt, knorm = kpoints.get_dirs(HS_dir, poscar, k_unit)
    if norm == None :
        if len(sp) == 0:
            norm=mpl.colors.Normalize(0,1)
        else :
            norm=mpl.colors.Normalize(-1,1)
    if ax == None:
        ax = plt.subplot()
    BANDS = Band_Structure_interp(path,Elim = Elim,norm=norm, Ef = Ef_SCF, orb=orb, atom=at, spin = sp, 
                                  interp_kind=interp_kind, HS_dir = HS_dir,cmap = cmap, BANDS = BANDS, picker = picker, Disp_BANDS = Disp_BANDS, 
                                  kpoints = kpoints, kpt = kpt, knorm =knorm, k_unit = k_unit); 
    
    if plot_ticks :
        ticks, label = search_sym_lines_and_plot(ax, kpoints, kpt, knorm = knorm, k_unit = k_unit)
        plt.hlines(0, *ax.get_xlim(), color = 'black', linestyle = ':')
    return BANDS, ax
    
def routine_BS_d(path, Elim, Ef = None, HS_dir = None, interp_kind = None,  seuil = 0.2, name = None, legend = True, ax = None, kpt = None, poscar = None, fig = None, new_fig = False):
    """
    (Vielle version, il vaut mieux utiliser celle défini dans lib_display_BS)
    
    
    Permet de tracer la stucture de bande et faire des projections sur les orbitales d et p.
    Actuellement le programme trace en fonction de l'index des kpoints et non la norme
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi (si None Ef est déterminé par le OUTCAR ou le DOSCAR du path)
    - HS_dir : pas encore implémenté
    - interp_kind : type d'interpolation si None --> pas d'interpolation ('cubic' / 'linear' / None)
    - seuil : seuil pour les projections.
    - name : Nom de la figure.
    - Legend = True : affiche la légende.
    - ax = None : axe matplotlib utilisé.
    
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        routine_BS_d(path, Elim = [Emin, Emax])
    """
    if name == None: 
        name = 'BS orb decomp'
    if fig == None :
        fig = plt.figure(name)
        new_fig = True
    norm = Normalize(seuil-0.1,seuil)
    if Ef == None :
        Ef = get_Ef(path, fileformat = 'PROCAR')
    ORBS = {'dxy' : cmap_r, 'dyz' : cmap_g, 'dxz' : cmap_p, 'dz2' : cmap_b, 'x2-y2' : cmap_c, 'p' : cmap_o}
    kpoints = KPOINTS(path, fileformat = 'PROCAR')
    if type(kpt) == type(None) :
        kpt, knorm = kpoints.remove_double(kpoints.data.index, np.zeros(len(kpoints)))
    Disp_BANDS = get_bands_limit(path, Elim, Ef)
    
    BANDS, ax = axes_BS_interp(path,Elim, Ef=Ef, orb = 'tot', HS_dir = HS_dir, interp_kind=interp_kind, cmap = cmap_n, norm = norm, ax = ax, figname=name, kpoints = kpoints, kpt = kpt, Disp_BANDS = Disp_BANDS, poscar = poscar, fig = fig, new_fig = new_fig)
    custom_lines = [Line2D([0], [0], color=cmap_n(0), label = 'tot')]
    for orb in ORBS.keys():
        cmap = ORBS[orb]
        custom_lines.append(Line2D([0], [0], color=cmap(0), label = orb))
        if orb == 'p':
            orb = 'px+py+pz'
        axes_BS_interp(path,Elim, Ef=Ef, orb = orb, HS_dir = HS_dir, interp_kind=interp_kind, cmap = cmap, norm = norm, ax = ax, BANDS = BANDS, figname=name, plot_ticks=False, new_fig = False, picker = None, kpoints = kpoints, kpt = kpt,Disp_BANDS = Disp_BANDS, poscar = poscar, fig = fig)
    if legend :
        if new_fig == True :
            #fig.legend(handles=custom_lines, bbox_to_anchor=(0,1.0,1,0.2), ncol = len(ORBS.keys())+1, mode = 'expand')   
            fig.legend(handles=custom_lines, loc = 'upper center', ncol = len(ORBS.keys())+1, mode = 'expand')   
    return fig
    
    
def get_bands_from_wannier(path, spin = '') :
    """
    Permet d'obtenir les bandes calculé en Wannier
    Renvoi un tableau contenant les différentes bandes et un tableau contenant la norme des k-points
    - path : chemin vers le dossier où wannier90_band.dat est present
    - spin : '' / 'up' / 'dn' si calcule en spin
    Utilisation :
    
    path = "chemin vers le dossier où wannier90_band.dat est present"
    BANDS_w, knorm = get_bands_from_wannier(path)
    
    # Plot :
    Ef = _
    Elim = [E_min, E_max]
    figname = "PLOT Wannier"
    fig = plt.figure(figname)
    
    for i in BANDS_w :
        plt.plot(knorm, i - Ef, color = 'blue')
        
    # pour avoir les points de HS :
    kpoints = KPOINTS(path, fileformat = 'wannier')
    ticks, labels = kpoints.ticks_for_plot(np.arange(len(knorm)), knorm)
    plt.xticks(knorm[ticks], labels)
    plt.vlines(knorm[ticks], *Elim, color = 'black')
    
    plt.ylim(*Elim)
    plt.xlim(knorm[0], knorm[-1])
    
    
    """
    BANDS = []
    KPTS = []
    BANDS.append([])
    if 'up' in spin :
        spin = '.up'
    elif 'dn' in spin :
        spin = '.dn'
    b = 0
    with open(path+'wannier90'+spin+'_band.dat', 'r') as f :

        for l in f:
            try :
                KPTS.append(float(l.split()[0]))
            except :
                print('kpts done')
                break
    with open(path+'wannier90'+spin+'_band.dat', 'r') as f :
        
        for l in f:
            try :
                BANDS[b].append(float(l.split()[1]))
            except :
                b +=1
                BANDS.append([])
    BANDS.remove([])
    
    return np.array(BANDS), np.array(KPTS)
      
#def routine(path, Ef = None,at = 'tot', orb = 'tot', sp = '', cmap = plt.get_cmap('seismic'), path_Ef = None, HS_dir = None, interp_kind = 'cubic' , ax = None, Elim = [-0.5,0.5], norm = None):
#fig = plt.figure()
#fig.canvas.mpl_connect('pick_event', onpick)       
         
class Bands():
    class BAND():
        def __init__(self, energy, occ) :
            self.bands = energy
            self.occ = occ
            
        def __repr__(self):
            string = "%5s   %8s   %5s\n" %("kpt", "Energy", "occ.")
            for kpt in range(len(self.bands)):
                string += "%5d   % 4.4f   %1.3f\n" %(kpt, self.bands[kpt], self.occ[kpt])
            return string
    
    def __init__(self, path, spin = None):
        self.ISPIN = get_ISPIN(path)
        self.path = path
        self.nbands = int(grep('NBANDS', path+'OUTCAR')[0][0].split()[-1])
        if self.ISPIN == 1 :
            self.spin = ''
        elif spin == None :
            self.spin = '_up'
        else :
            self.spin = spin
        self.load_outcar()
    
    def load_outcar(self):
        bands_info = np.array(get_band_info(self.path,spin = self.spin))
        self.bands = bands_info
    
    def __getitem__(self, b):
        return self.bands[:,b]
    
    def get_bands_in_Elim(self, Elim, dE = 0):
        bands_account = np.array([])
        for b in range(self.nbands) : 
            band = self.bands[0][b]
            if np.any(band >= Elim[0]+dE) and np.any(band <= Elim[1]+dE):
                bands_account = np.append(bands_account,b)
        return np.int0(bands_account)

    def plot(self, Elim, Ef, kpt = None, knorm = None, bands_account = None,
             ax = None, 
             **kwargs):
        """
        """
        if type(bands_account) == type(None):
            bands_account = self.get_bands_in_Elim(Elim, dE = Ef)
        if type(kpt) == type(None) :
            kpt = np.arange(len(self.bands[0][0]))
        if type(knorm) == type(None) :
            knorm = kpt
        if ax == None :
            ax = plt.gca()
            
        for b in bands_account :
            #print (kpt, knorm,self.bands[0][b][kpt])
            ax.plot(knorm, self.bands[0][b][kpt]-Ef, label = str(b), **kwargs)
        return ax
            
    def routine(self, Elim, Ef = None, HS_lines = 'auto', ax = None,**kwargs):
        ref = time.time()
        try :
            self.kpt
            self.knorm
            self.kpoints = KPOINTS(self.path)
        except :
            self.kpoints = KPOINTS(self.path)
            self.poscar = POSCAR(self.path)
            kpoints = self.kpoints
            poscar = self.poscar
                
            fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
            if type(HS_lines) == str :
                if HS_lines == 'auto' :
                    HS_lines = list(kpoints.auto_lines())
            self.kpt, self.knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz)
        print(time.time()-ref)
        if ax == None :
            ax = plt.gca()
        self.plot(Elim, Ef,self. kpt, self.knorm, **kwargs)
        ticks, label = new_ticks(self.kpoints,self.kpt,self.knorm)
        #ax.xticks(knorm[np.int0(ticks)],label)
        ax.set_xticks(self.knorm[np.int0(ticks)])
        ax.set_xticklabels(label)
        ax.vlines(self.knorm[np.int0(ticks)],*Elim, label = 'TICKS')
        ax.set_xlim(self.knorm[0], self.knorm[-1])
        ax.set_ylim(*Elim)
        ax.hlines(0,*ax.get_xlim(), color = 'black', linestyle = ':')
        
            
#path = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/constrain/a388_c416/SOC/BS/"            
path = "/home/gosteau/Documents/These/2eme année/STO surface/sym/centered/1.5/BS/"           
            
            
            
            
