#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 10:07:12 2018

@author: gosteau
"""

from lib_gen import *
from lib_poscar import *
ORB_F = True
# A immplementer :
"""
1) recuperér le DOSCAR par atome et le stocker (ecriture fichier) self.doscar (pd.DataFrame)
 i) dépend de LORBIT (s, p,d ... ou s, px,py,...)       --> self.lorbit
 ii) dépend des orbitales des atomes : orbitales f ...  --> self.orbitals
 iii) dépend du spin                                    --> self.spin
 iv) atome en question                                  --> self.atom (pd.DataFrame)

"""
ini_dos = 6
s = ['s']
p = ['py', 'pz', 'px']
d = ['dxy' , 'dyz' , 'dz2' , 'dxz' , 'dx2-y2']
f = ['f1' , 'f2' , 'f3', 'f4', 'f5' , 'f6' , 'f7']
ORBS = {'s' : s,'p':p,'d':d,'f':f}
#default_orb = 'yes'

def get_ISPIN(path):
    return int(grep("ISPIN",path + "OUTCAR")[0][0].split()[2])

def get_LORBIT(path):
    return int(grep("LORBIT", path+"OUTCAR")[0][0].split()[2])

def get_nDOS(path):
    return int(grep("NEDOS", path+"OUTCAR")[0][0].split("=")[1].split()[0])
    
def get_LSORBIT(path):
    LSORBIT = grep("LSORBIT", path+"OUTCAR")[0][0].split("=")[1].split()[0]
    if LSORBIT == 'T':
        return True
    else :
        return False
    
def get_orbitals(path) :
    def l_to_orb(l):
        if l == 0:
            return 's'
        elif l == 1:
            return 'p'
        elif l == 2:
            return 'd'
        elif l == 3 and ORB_F == True:
            return 'f'
        else :
            return 'none'
        
    l_pos = grep('Description',path+'OUTCAR')[1]
    n_pos = np.array(grep('local pseudopotential read in',path+'OUTCAR')[1])-np.array(l_pos)
    L = [[],[]]
    
    for i in range(0,len(l_pos)):
        ligne = l_pos[i]
        n = n_pos[i]
        data = pd.read_csv(path+'OUTCAR',skiprows = ligne,nrows=n-2, delim_whitespace=True)
        for l in data['l'].sort_values() :
            if l not in L[0] :
                
                tmp = l_to_orb(l)
                if tmp != 'none' :
                    L[0].append(l)
                    L[1].append(tmp)
    return L[1]



class DOSCAR: 
    def default_orbs(self):
        LORBIT = self.lorbit
        if LORBIT == 11 :
            default_orb = []
            for l in self.outcar_orbitals:
                default_orb += ORBS[l]
            for l in self.outcar_orbitals:
                if l != 's' :
                    default_orb += [l]
            default_orb += ['tot']
        else :
            default_orb = self.outcar_orbitals+['tot']
        return default_orb

    def head(self,path, atom, orbitals) :
        spin=self.spin
        HEAD=['Energy']
        if spin == 2 :
            spin_st = ['_up' , '_dn' ]
        else :
            spin_st = [''] 
        if self.lsorbit :
            spin_st = ['_stot','_sx','_sy','_sz']
        for l in orbitals :
                for sp in spin_st:
                    HEAD.append(str(atom)+'_'+l+sp)
        return HEAD
    
    # DOS_plan head
    def data_head(self,path,atoms, orbitals) :
        
        spin=self.spin
        HEAD=['Energy']
        if spin == 2 :
            spin_st = ['_up' , '_dn' ]
        else :
            spin_st = [''] 
        if self.lsorbit :
            spin_st = ['_stot','_sx','_sy','_sz']
        for atom in atoms :
            HEAD = HEAD + self.head(path,atom,orbitals)[1:]
        return HEAD    
    
    def get_energy_range(self,path,Elim= None):
        energy = pd.read_csv(path+'DOSCAR',names=['Energy'],usecols=['Energy'],skiprows=6,nrows=self.n_dos,delim_whitespace=True)
        
        if Elim != None :
            energy = energy.loc[energy['Energy']-self.Ef >= Elim[0]]
            energy = energy.loc[energy['Energy']-self.Ef <= Elim[1]]
        return energy
        
    
    def get_at_DOS(self,path, atom, get_data = True):
        timeref = time.time()
        filename = path+"Atoms_dos/"+"ATOM_"+str(atom)
        if not os.path.isdir(path+"Atoms_dos") :
            os.mkdir(path+'Atoms_dos')
        if not os.path.isfile(filename):
            HEAD = self.head(path,atom,self.default_orbitals)
            data = pd.read_csv(path+'DOSCAR',names = HEAD, skiprows=ini_dos+(self.n_dos+1)*atom, nrows=self.n_dos,delim_whitespace = True)
            spin = self.spin
            if spin == 2 :
                spin_st = ['_up' , '_dn' ]
            else :
                spin_st = [''] 
            if self.lsorbit :
                spin_st = ['_stot','_sx','_sy','_sz']
                
            LORBIT = self.lorbit
            for sp in spin_st :
                data[str(atom)+'_'+'tot'+sp] = 0.
                if LORBIT == 11 :
                    for l in self.outcar_orbitals :
                        if l != 's':
                            data[str(atom)+'_'+l+sp] = 0.
                            for m in ORBS[l]:
                                if np.isnan(data[str(atom)+'_'+m+sp][0]) :
                                    data[str(atom)+'_'+m+sp] = 0.
                                data[str(atom)+'_'+l+sp] += data[str(atom)+'_'+m+sp]
                                data[str(atom)+'_'+'tot'+sp] += data[str(atom)+'_'+m+sp]
                        else :
                            data[str(atom)+'_'+'tot'+sp] += data[str(atom)+'_'+'s'+sp]
                else :
                    for l in self.outcar_orbitals:
                        data[str(atom)+'_'+'tot'+sp] += data[str(atom)+'_'+l+sp]
                    
            data.set_index('Energy').to_csv(filename, sep=' ')
            
            #pkl.dump(data,file(filename,'wb'))
        if get_data :
            HEAD = self.head(path,atom,self.used_orbitals)
            
            if self.rest_energy :
                datahead = get_row(1,filename).split()
                index_inf = self.energy.index[0]+1
                nrows = len(self.energy)
                #index_inf = self.energy.index[0]
                #nrows = len(self.energy)-1
                data = pd.read_csv(filename,names = datahead,usecols=HEAD,skiprows=index_inf,nrows=nrows, delim_whitespace = True)
                #data = pkl.load(file(filename,'rb'))[HEAD].loc[index_inf:nrows]
            else :
                data = pd.read_csv(filename,usecols=HEAD, delim_whitespace = True)
                #data = pkl.load(file(filename,'rb'))[HEAD]
            
            #print "Atom :" +str(atom)+' : '+ str(time.time()-timeref)
            return data
        
    def get_atoms_DOS(self,path, atoms):
        PLANE_DOS = pd.DataFrame(columns = self.data_head(path,atoms,self.used_orbitals))
        if self.spin == 2 :
            spin_st = ['_up' , '_dn' ]
        else :
            spin_st = [''] 
        if self.lsorbit :
            spin_st = ['_stot','_sx','_sy','_sz']
        for p in atoms :
            data = self.get_at_DOS(path, p)
            for c in data.columns[1:] :
                PLANE_DOS[c] = data[c]
            
            for orb in self.used_orbitals :
                #for orb in princ_orbs.keys()+['tot'] :
                for sp in spin_st :
                    if p == atoms[0] :
                        PLANE_DOS[orb+sp] = 0.
                    PLANE_DOS[orb+sp] += PLANE_DOS[str(p)+'_'+orb+sp]
             
        PLANE_DOS['Energy'] = data['Energy']        
        return PLANE_DOS
    
    def get_DOS(self, path):
        if self.spin == 1 :
            columns = ['Energy', 'dos', 'sum']
        else :
            columns = ['Energy', 'dos_up', 'dos_dn','sum_up', 'sum_dn']
        data = pd.read_csv(path+"DOSCAR",names=columns, delim_whitespace = True,skiprows = 6, nrows = self.n_dos)
        return data
            
    def SbH(self, orbital, seuil = 0.1):
        try :
            #print np.array(self.data.loc[self.data['Energy']-self.Ef > 0]['Energy'] -self.Ef)
            #print np.array(self.data.loc[self.data['Energy']-self.Ef > 0].loc[self.data[orbital] < seuil]['Energy'] -self.Ef)
            Ec = np.array(self.data.loc[self.data['Energy']-self.Ef > 0].loc[self.data[orbital] < seuil]['Energy'] -self.Ef)[-1]
            
            #Ec = np.array(self.data.loc[self.data[orbital]-self.Ef > 0].loc[self.data[orbital] < seuil]['Energy'] -self.Ef)[-1]
        except :
            Ec = np.nan
        try :
            #Ev = np.array(self.data.loc[self.data['Energy']-self.Ef < 0].loc[self.data[orbital] < seuil]['Energy'] -self.Ef)[0]
            Ev = np.array(self.data.loc[self.data['Energy']-self.Ef < 0].loc[self.data[orbital] < seuil]['Energy'] -self.Ef)[0]
        except :
            Ev = np.nan
        return Ec,Ev
        
    
    def __init__(self,path, atoms =-1, orbitals = None, Ef = None,Elim = None,energy_range = None, LORBIT = None):
        if LORBIT == None :
            LORBIT = get_LORBIT(path)
        SPIN = get_ISPIN(path)
        n_dos = get_nDOS(path)
        if Ef == None :
            Ef = get_Ef(path,same_dir=True, fileformat = 'DOSCAR')
        self.lorbit = LORBIT
        self.spin = SPIN
        self.Ef = Ef
        self.n_dos = n_dos
        self.lsorbit = get_LSORBIT(path)
        self.outcar_orbitals = get_orbitals(path)
        self.default_orbitals = self.default_orbs()
        if Elim != None or type(energy_range) != type(None) :
            self.rest_energy = True
        else :
            self.rest_energy = False
            
        if orbitals == None :
            self.used_orbitals = self.default_orbitals
        else :
            self.used_orbitals = orbitals
        
        if type(energy_range) != type(None):
            self.energy = energy_range
        else :
            self.energy = self.get_energy_range(path,Elim)
        if type(atoms) == int:
            if atoms == -1 :
                self.data = self.get_DOS(path)
            else :
                self.data = self.get_at_DOS(path,atoms)
                self.atoms = [atoms]
        else :
            self.data = self.get_atoms_DOS(path,atoms)
            self.atoms = atoms
            
        

class CHGDENS2:
    def get_CHGDENS(self,planar_DOSCAR, PLANES, sign=1):
        orbs = []
        for orb in self.orbitals :
            for sp in self.uspin :
                orbs.append(orb+sp)
        Data = pd.DataFrame(index = range(0,len(PLANES)),columns = orbs)
        if sign == 1 :
            liminf = self.Ef-self.dE
            limsup = self.Ef
        else :
            liminf = self.Ef
            limsup = self.Ef+self.dE
        for p in range(0,len(PLANES)) :
            DOS = planar_DOSCAR[p]
            DOS = DOS.loc[DOS['Energy']> liminf]
            DOS = DOS.loc[DOS['Energy']<= limsup]
            for orb in self.orbitals :
                for sp in self.uspin :
                    Data[orb+sp][p] = sign*scp.simps(DOS[orb+sp], DOS['Energy'])
        return Data
    
    def get_CHGDENS_LAOSTO(self):
        poscar = POSCAR(self.path)
        PLANE = poscar.get_plane()
        n_at = np.array(poscar.data.set_index('At').loc['Ti'].reset_index()['n_at'])[0]
        IF1 = get_IF(n_at,PLANE)
        n_at = np.array(poscar.data.set_index('At').loc['Ti'].reset_index()['n_at'])[-1]
        IF2 = get_IF(n_at,PLANE)
        orbs = []
        for orb in self.orbitals :
            for sp in self.uspin :
                orbs.append(orb+sp)
        CHGDENS = pd.DataFrame(columns = orbs)
        DOSCAR_planar = [[]]*len(PLANE)
        for p in range(0,len(PLANE)):
            DOSCAR_planar[p] = DOSCAR(self.path,PLANE[p],orbitals=self.orbitals,Ef= self.Ef).data
        tab = [[0,IF1],[IF1,IF2+1],[IF2+1,-1]]
        sign = 1
        for i in tab :  
            sign *= -1
            if i[1] != -1 :
                tmp = self.get_CHGDENS(DOSCAR_planar[i[0]:i[1]], PLANE[i[0]:i[1]],sign=sign)
            else : 
                tmp = self.get_CHGDENS(DOSCAR_planar[i[0]:], PLANE[i[0]:],sign=sign)
            CHGDENS = CHGDENS.append(tmp, ignore_index = True)
        return CHGDENS
    
    def __init__(self, path, Ef = None, dE = 0.5, orbitals = ['tot'],uspin = ['']):
        self.dE = dE
        self.path = path
        if Ef == None :
            Ef = get_Ef(path,same_dir=True)
        SPIN = get_ISPIN(path)
        self.Ef = Ef
        self.orbitals = orbitals
        self.spin = SPIN
        if self.spin == 1:
            self.uspin = ['']
        else :
            self.uspin = []
            count = 0
            for s in uspin :
                if s in ['_up','_dn']:
                    count += 1
                    self.uspin += [s]
            if count == 0 :
                self.uspin = ['_up']
        
        self.data = self.get_CHGDENS_LAOSTO()

class CHGDENS:
    # NEW !!!!!!!!   
        
    def get_CHGDENS(self, atoms, sign=1, doscar = None):
        def get_power_2(x):
            i=0
            while 2**i < x:
                i+=1
            dx=x-2**(i-1)
            return dx
        orbs = []
        for orb in self.orbitals :
            for sp in self.uspin :
                orbs.append(orb+sp)
        
        #doscar = DOSCAR(self.path,atoms,orbs,Elim = [-self.dE,self.dE])
        
        if sign == 1:
            echo = True
            if doscar == None :
                doscar = DOSCAR(self.path,atoms,self.orbitals,Elim = [-self.dE,0])
            else :
                a= doscar.data['Energy']-doscar.Ef 
                b = a.loc[a >= -self.dE]
                c = b.loc[b <= 0]
                doscar.data = doscar.data.loc[c.index]
            if len(doscar.data)%2 == 0:
                doscar.data = doscar.data.loc[19:]
        else :
            echo = False
            if doscar == None :
                doscar = DOSCAR(self.path,atoms,self.orbitals,Elim = [0,self.dE])
            else :
                a= doscar.data['Energy']-doscar.Ef 
                b = a.loc[a >= 0]
                c = b.loc[b <= self.dE]
                doscar.data = doscar.data.loc[c.index]
            if len(doscar.data)%2 == 0:
                doscar.data = doscar.data.loc[:len(doscar.data)-20]
                
        if self.auto :
            a= doscar.data['tot'+self.uspin[0]] >= self.zeros
            doscar.data = doscar.data.loc[a]
            #print doscar.data[['Energy','tot']]
            l = len(doscar.data)
            
            dx = get_power_2(l) -1
            if self.methods == 'romb':
                if sign == 1:
                    doscar.data = doscar.data.reset_index().loc[dx:]
                else :
                    doscar.data = doscar.data.reset_index().loc[:l-dx-1]
                l2 = len(doscar.data)
        Data = {}
        for orb in orbs :
            if l == 0 :
                Data[orb] = 0
            else :
                if self.methods == 'simps' :
                    Data[orb] = sign*scp.simps(np.array(doscar.data[orb]), np.array(doscar.data['Energy']))
                elif self.methods == 'trapz' :
                    Data[orb] = sign*scp.trapz(np.array(doscar.data[orb]), np.array(doscar.data['Energy']))
                elif self.methods == 'romb' :
                    Data[orb] = sign*scp.romb(np.array(doscar.data[orb]), np.array(doscar.data['Energy']))
                elif self.methods == 'cumtrapz' :
                    Data[orb] = sign*scp.cumtrapz(np.array(doscar.data[orb]), np.array(doscar.data['Energy']))
        #Data= {orb : sign*scp.simps(np.array(doscar.data[orb]), np.array(doscar.data['Energy'])) for orb in orbs}
        #print Data
        
        #print atoms, sign, Data['tot'],echo
        return Data
    
    
    def get_CHGDENS_LAOSTO(self,PLANE,planar_dos = None):
        if planar_dos != None :
            if len(PLANE) != len(planar_dos) :
                message = "lenght of PLANE and planar_dos doesn't correspond :\n"+"len(PLANE) = "+str(len(PLANE))+"\n"+"len(planar_dos) = "+str(len(planar_dos))
                ERROR(message)
            
            for i in range(0,len(PLANE)) :
                for at in PLANE[i] :
                    if at not in planar_dos[i].atoms :
                        message = str(at)+' from ' +str(i)+' PLANE not in '+str(i)+' planar_dos\nPLANE = '+str(PLANE[i])+'\n'+'planar_dos.atoms = '+str(planar_dos[i].atoms)
                        ERROR(message)
        
        orbs = []
        for orb in self.orbitals :
            for sp in self.uspin :
                #print orb+sp
                orbs.append(orb+sp)
                if planar_dos != None:
                    if orb+sp not in planar_dos[0].used_orbitals :
                        message = "input orbitals don't correspond to orbitals of planar_dos"
                        ERROR(message)
        Data = pd.DataFrame(columns = orbs)
        
        
        
        for i in range(0,len(PLANE)) :
            if i in self.holes :
                sign = -1
            else :
                sign = 1
            if planar_dos != None :
                tmp = self.get_CHGDENS(PLANE[i],sign,planar_dos[i])
            else :
                tmp = self.get_CHGDENS(PLANE[i],sign)
            Data = Data.append(tmp,ignore_index=True)
        return Data
        
    def __init__(self, path, PLANE, holes = [], Ef = None, dE = None, dE_max = 0.5, orbitals = ['tot'],uspin = [''], zeros = 1e-3, methods = 'simps', planar_dos = None, SUM_spin= True):
        if dE == None :
            self.auto = True
            dE = dE_max
        else :
            self.auto = False
        self.dE = dE
        self.methods = methods
        self.holes = holes
        self.zeros = zeros
        self.path = path
        if Ef == None :
            Ef = get_Ef(path,same_dir=True)
        SPIN = get_ISPIN(path)
        self.Ef = Ef
        self.orbitals = orbitals
        self.spin = get_ISPIN(path)
        self.SUM_spin = SUM_spin
        self.lsorbit = get_LSORBIT(path)
        if self.spin == 1:
            self.uspin = ['']
        elif self.spin == 2 :
            self.uspin = []
            count = 0
            for s in uspin :
                if s in ['_up','_dn']:
                    count += 1
                    self.uspin += [s]
            if count == 0 :
                self.uspin = ['_up']
        if self.spin == 2 and self.SUM_spin :
            self.uspin = ['_up','_dn']
        if self.lsorbit :
            self.uspin = []
            count = 0
            for s in uspin :
                if s in ['_stot','_sx','_sy','_sz']:
                    count += 1
                    self.uspin += [s]
            if count == 0 :
                self.uspin = ['_stot']
        #print self.uspin
        
        self.data = self.get_CHGDENS_LAOSTO(PLANE,planar_dos)
        if self.SUM_spin and self.spin == 2 :
            for orb in self.orbitals :
                self.data[orb] = self.data[orb+'_up']+ self.data[orb+'_dn']



def chgdens_display(chgdens, name = 'CHGDENS',figsize = (8,6),orbs = {'tot':'blue','dxy':'red','dxz':'purple','dyz':'green'},
uspin = [''],range_plane = None,label = [], vlines = {},colorsbox= None,surf = None, methods = 'simps',ylim = None):
    if range_plane == None :
        range_plane = range(0,len(chgdens.data))
    majorLocator = MultipleLocator(1)
    minorLocator = MultipleLocator(2)
    fig = plt.figure(name,figsize)
    ax = plt.subplot()
    uspin = chgdens.uspin
    print(chgdens.data.columns)
    if chgdens.SUM_spin :
        uspin = ['']
    for orb in orbs.keys() :
        for sp in uspin :
            print(orb+sp)
            if orb+sp in chgdens.data.columns :
                plt.plot(chgdens.data.index[range_plane]+1,chgdens.data[orb+sp][range_plane], marker = '+', linestyle = '-', color = orbs[orb])
            elif len(orb.split('+')) > 1 :
                chgdens.data[orb] = 0
                for l in orb.split('+') :
                    chgdens.data[orb+sp] += chgdens.data[l+sp]
                plt.plot(chgdens.data.index[range_plane]+1,chgdens.data[orb+sp][range_plane], marker = '+', linestyle = '-', color = orbs[orb])
    if ylim == None :
        ylim = plt.ylim()
    else :
        plt.ylim(*ylim)
    xlim = (range_plane[0]+1,range_plane[-1]+1)
    plt.xlim(*xlim)
    for v in vlines.keys() :
        plt.vlines(v+1,*ylim,colors='black',linestyles=vlines[v])
    plt.hlines(0,*xlim,colors='black', linestyles='--')
    
    if colorsbox != None :
        for i in range(0,len(colorsbox['box'])):
            x0 = colorsbox['box'][i][0]+1
            dx = colorsbox['box'][i][1]+1-x0
            y0 = ylim[0]
            dy = ylim[1]-ylim[0]
            color = colorsbox['color'][i]
            alpha = colorsbox['alpha']
            rec =plt.Rectangle((x0,y0),dx,dy,color=color,alpha=alpha)
            plt.gca().add_patch(rec)
    
    tab = range(range_plane[0]+1,range_plane[-1]+2,2)
    plt.xticks(tab,label)
    ax.xaxis.set_minor_locator(minorLocator)
    ax.tick_params(which='both', width=1)
    ax.tick_params(which='major', length=7)
    ax.tick_params(which='minor', length=4)
    
    locs,label = plt.yticks()
    plt.ylim(*ylim)
    ax.set_yticks(locs,label)
    ax.set_ylim(*ylim)
    ax.set_ylabel('N(z) (states/f.u.)')
    if surf != None :
        ax2=plt.twinx()
        ax2.set_yticks(locs)
        ax2.set_ylim(*ylim)
        ax2.set_yticklabels(['%.1e' % (i/surf) for i in locs])
        ax2.set_ylabel("N(z) (states/cm2)", color = 'black')
    
    return fig




path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/GGA/StageM2 lattice/"
path2 = "/home/gosteau/Documents/TMP_remi/FeSTO/vac/"
path = "/home/gosteau/Documents/Stage M2/Résultats/20STO/0.000001/RUN2_EFIELD_0.000001/"
path_dos = "/home/gosteau/Documents/Stage M2/Résultats/20STO/0.0024/RUN_EFIELD_0.0024/"



def routine_chgdens(path,zeros = 2e-2,surf = True, name = 'CHGDENS', figsize = (10,6), orbitals =  None,orbs = None, SUM_spin = True,methods = 'simps', ylim=None, poscar = None) :
    if poscar == None :
        poscar = POSCAR(path);
    lorbit = get_LORBIT(path)
    if orbitals == None and orbs == None :
        if lorbit == 10 :
            orbitals = ['tot','d','p']
            orbs={'tot':'blue', 'd':'green', 'p':'red'}
        elif lorbit == 11 :
            orbitals = ['tot','dxy','dxz','dyz','px','py','pz']
            orbs={'tot':'blue', 'dxy+dyz+dxz':'green', 'px+py+pz':'red'}
            orbs={'tot':'blue', 'dyz+dxz':'green', 'dxy':'red'}
        else :
            orbitals = ['tot']
            orbs={'tot':'blue'}
            
    plane = poscar.get_plane(); 
    n_at = poscar.data.sort_values('z').set_index('At').loc['Ti']['n_at'][0]
    IF1 = get_IF(n_at,plane)
    n_at = poscar.data.sort_values('z').set_index('At').loc['Ti']['n_at'][-1]
    IF2 = get_IF(n_at,plane)
    holes = range(0,IF1)+range(IF2+1,len(plane)); 
    #holes = []
    chgdens = CHGDENS(path,plane,holes, orbitals = orbitals,zeros = zeros, SUM_spin = SUM_spin,methods= methods)
    range_plane = range(0,IF1)+range(IF1,IF2+1,2)+range(IF2,len(chgdens.data))
    #range_plane = range(IF1,IF2+1,2)
    label = []
    for i in range(0,len(chgdens.data),2):
    #for i in range_plane :
        if i <= (IF1+IF2)/2 :
            label += [(IF1-i)/2]
        else : 
            label += [(i-IF2)/2]
    
    vlines = {IF1 : '-',IF2 : '-', (IF1+IF2)/2 : '--'}
    colorsbox = {'box' : [[0,IF1],[IF1,IF2],[IF2,len(chgdens.data)]], 'color' : ['purple','green','purple'], 'alpha' : 0.2}
    if surf :
        surf = poscar.diag_lattice[0]*poscar.diag_lattice[1]*(1e-8)**2
    else :
        surf = None
    return chgdens, chgdens_display(chgdens,name = name,figsize = figsize,range_plane = range_plane, vlines=vlines, label = label, colorsbox=colorsbox, surf = surf, orbs = orbs,ylim = ylim)

def routine_chgdens_asym(path,zeros = 2e-2,surf = True, name = 'CHGDENS', figsize = (10,6), orbitals =  None,orbs = None, SUM_spin = True,methods = 'simps', ylim=None) :
    poscar = POSCAR(path);
    lorbit = get_LORBIT(path)
    if orbitals == None and orbs == None :
        if lorbit == 10 :
            orbitals = ['tot','d','p']
            orbs={'tot':'blue', 'd':'green', 'p':'red'}
        elif lorbit == 11 :
            orbitals = ['tot','dxy','dxz','dyz','px','py','pz']
            orbs={'tot':'blue', 'dxy+dyz+dxz':'green', 'px+py+pz':'red'}
            orbs={'tot':'blue', 'dyz+dxz':'green', 'dxy':'red'}
        else :
            orbitals = ['tot']
            orbs={'tot':'blue'}
            
    plane = poscar.get_plane(); 
    n_at = poscar.data.sort_values('z').set_index('At').loc['Ti']['n_at'][0]
    IF1 = get_IF(n_at,plane)
    n_at = poscar.data.sort_values('z').set_index('At').loc['Ti']['n_at'][-1]
    IF2 = get_IF(n_at,plane)
    n_at = poscar.data.sort_values('z').set_index('At').loc['Sr']['n_at'][-1]
    IF2 = get_IF(n_at,plane)
    holes = range(0,IF1)+range(IF2+1,len(plane)); 
    #holes = []
    chgdens = CHGDENS(path,plane,holes, orbitals = orbitals,zeros = zeros, SUM_spin = SUM_spin,methods= methods)
    range_plane = range(0,IF1)+range(IF1,IF2+1,2)+range(IF2,len(chgdens.data))
    #range_plane = range(IF1,IF2+1,2)
    label = []
    for i in range(0,len(chgdens.data),2):
    #for i in range_plane :
        label += [(IF1-i)/2]
    n_at = poscar.data.sort_values('z').set_index('At').loc['Sr']['n_at'][-1]
    IF2 = get_IF(n_at,plane)
    vlines = {IF1 : '-',IF2 : '-'}
    colorsbox = {'box' : [[0,IF1],[IF1,IF2],[IF2,len(chgdens.data)]], 'color' : ['purple','green','purple'], 'alpha' : 0.2}
    if surf :
        surf = poscar.diag_lattice[0]*poscar.diag_lattice[1]*(1e-8)**2
    else :
        surf = None
    return chgdens, chgdens_display(chgdens,name = name,figsize = figsize,range_plane = range_plane, vlines=vlines, label = label, colorsbox=colorsbox, surf = surf, orbs = orbs,ylim = ylim)



path1 = "/home/gosteau/Documents/Stage M2/Résultats/20STO/wo field/DOS_more_kpt/"
path2 = "/home/gosteau/Documents/Stage M2/Résultats/20STO/SOC/vanilla/DOSCAR_SOC_SCF_3/"
path = "/home/gosteau/Documents/Stage M2/Résultats/20STO/0.0018/RUN_EFIELD/"
path = "/home/gosteau/Documents/Stage M2/Résultats/20STO/NELECT/1101.74/"
#path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/sym/vanilla/test/GGA_COMPAT/+U/dos_res_2/"
#path='/home/gosteau/Documents/TMP_remi/LAOSTO_Ir/without_Ir/'



path = "/home/gosteau/Documents/TMP_remi/LAOSTO_Ir/"
IF = ['IF','IF-2','IF+1','IF+3','S']

IF = ['IF-3']

path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.8661/Relax/"
path2 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.90555/Relax/"
path3 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.98445/Relax/"
path4 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/4.0239/Relax/"
path_ref = '/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/sym/vanilla/test/GGA_COMPAT/+U/dos_res_2/'

path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/GGA/test_12STO/Ti/"
path2 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/GGA/test_12STO/Ti_pv/"
path3 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/GGA/test_12STO/LSDA/"


path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.8661/DOS/"
path1 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.90555/DOS/"
path2 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.98445/DOS/"
path3 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/4.0239/DOS/"

path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.8661/ISPIN2/DOS/"
path1 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.90555/ISPIN2/DOS/"
path2 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.98445/ISPIN2/DOS/"
path3 = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/4.0239/ISPIN2/DOS/"
path_ref = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.945/ISPIN2/DOS/"
methods = 'trapz'
ylim = None
ylim = [-0.07,0.03]
zeros = 1e-6

def plot_chgdens(chg,PATH, orbital):
    IF = 10
    Center = len(chgdens.data)/2
    dxy = []
    dxy_surf = []
    SURF = []
    for i in range(0,len(PATH)):
        C = chg[i]
        surf = PATH[i].split('/')[-3]
        print(i, surf)
        try :
            surf = float(surf)*1e-8
        except :
            surf = 3.945*1e-8
        if len(orbital.split('+')) != 1 :
            C.data[orbital] = 0
            for l in orbital.split('+')[1:]:
                C.data[orbital] += C.data[l]
        SURF.append(surf*1e8/3.945*100-100)
        dxy.append(np.sum(C.data[orbital][IF:Center:2]))
        dxy_surf.append(np.sum(C.data[orbital][IF:Center:2])/surf**2)
    fig1 = plt.figure(orbital)
    ax = plt.subplot()
    plt.plot(SURF,dxy, marker = '+', linestyle = '-')
    ax.set_ylabel('N(z) (states/f.u.)')
    ax.set_xlabel('strain (A)')
    fig2 = plt.figure(orbital+'norm')
    ax = plt.subplot()
    plt.plot(SURF,dxy_surf, marker = '+', linestyle = '-')
    ax.set_ylabel('N(z) (states/cm2)')
    ax.set_xlabel('strain (A)')
    return fig1, fig2

#chgdens,fig =routine_chgdens(path,zeros=zeros,SUM_spin=True,name='CHGDENS_'+path.split('/')[-4], methods = methods, ylim = ylim)
#chgdens1,fig1 =routine_chgdens(path1,zeros=zeros,SUM_spin=True,name='CHGDENS_'+path1.split('/')[-4], methods = methods, ylim = ylim)
#chgdens2,fig2 =routine_chgdens(path2,zeros=zeros,SUM_spin=True,name='CHGDENS_'+path2.split('/')[-4], methods = methods, ylim = ylim)
#chgdens3,fig3 =routine_chgdens(path3,zeros=zeros,SUM_spin=True,name='CHGDENS_'+path3.split('/')[-4], methods = methods, ylim = ylim)
#chgdens_ref,fig_ref =routine_chgdens(path_ref,zeros=zeros,SUM_spin=True,name='CHGDENS_'+path_ref.split('/')[-4], methods = methods, ylim = ylim)

#FIG = [fig,fig1,fig2,fig3,fig_ref]
"""
filename = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/Images/figures.pkl"
with bz2.BZ2File(filename,'rb') as f :    
    FIG = pkl.load(f)
"""
#chg = [chgdens,chgdens1,chgdens2,chgdens3,chgdens_ref]
filename = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/Images/chg.pkl"



def plot_chgdens(chg,PATH, orbital,ref = 2, x = None):
    if x == None :
        IF = 10
        Center = len(chg[ref].data)/2
    else :
        IF, Center = x
    dxy = []
    dxy_surf = []
    SURF = []
    
    C_ref = chg[ref]
    if len(orbital.split('+')) != 1 :
        C_ref.data[orbital] = 0
        for l in orbital.split('+')[1:]:
            C_ref.data[orbital] += C_ref.data[l]
    
    for i in range(0,len(PATH)):
        C = chg[i]
        surf = PATH[i].split('/')[-4]
        print(i, surf)
        try :
            surf = float(surf)*1e-8
        except :
            surf = 3.945*1e-8
        if len(orbital.split('+')) != 1 :
            C.data[orbital] = 0
            for l in orbital.split('+')[1:]:
                C.data[orbital] += C.data[l]
                
        RES = np.sum(C.data[orbital][IF:Center:2])/np.sum(C_ref.data[orbital][IF:Center:2])*100-100
        SURF.append(surf*1e8/3.945*100-100)
        dxy.append(RES)
        dxy_surf.append(np.sum(C.data[orbital][IF:Center:2])/surf**2)
    fig1 = plt.figure(orbital)
    ax = plt.subplot()
    plt.plot(SURF,dxy, marker = '+', linestyle = '-')
    xlim = plt.xlim()
    plt.xlim(*xlim)
    plt.hlines(0,*xlim,color = 'black', linestyle = ':')
    ax.set_ylabel('N(z) variation (%)')
    ax.set_xlabel('strain (A)')
    
    return fig1
"""
path = "/home/gosteau/Documents/Stage M2/Résultats/20STO/wo field/DOSCAR_20_STO/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/PBESOL/wo_U/DOS/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/20STO/constrain/3.945/Relax2/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/asym/PBESOL/8_STO/relax_higher/next/more_highter/DOS/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/verif_GGA_stage/GGA_COMPAT_T/exp/ISYM/Relax/DOS/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/verif_GGA_stage/GGA_COMPAT_T/exp/NELECT/test_1/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/6STO/PBESOL/Relax/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/verif_GGA_stage/GGA_COMPAT_T/exp/NELECT/relax_2/"
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/NP_IF/new_relax/"
name = '_IF_NP'
chgdens, fig = routine_chgdens(path, zeros = 1e-5, name = "CHGDENS"+name, methods = "trapz")
dxy = np.sum(chgdens.data.loc[10:50]['dxy'])/2
dmix = np.sum(chgdens.data.loc[10:50]['dxz']+chgdens.data.loc[8:50]['dyz'])/2
print dxy, dmix
"""
"""
path = "/home/gosteau/Documents/These/1ere année/LAOSTO/6STO/CHGDENS/"
name = '6STO'
chgdens, fig = routine_chgdens(path, zeros = 1e-5, name = "CHGDENS"+name, methods = "simps")

path = "/home/gosteau/Documents/These/1ere année/LAOSTO/7STO/PBESOL/D4h/CHGDENS/"
name = '7STO'
chgdens, fig = routine_chgdens(path, zeros = 1e-5, name = "CHGDENS"+name, methods = "simps")

"""

