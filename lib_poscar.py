#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 10:43:59 2018

@author: gosteau
"""

from lib_gen import *


LABEL_position = 'position of ions in fractional coordinates (direct lattice)'
LABEL_nion = 'ions per type'
LABEL_ion_type = 'VRHFIN'


class POSCAR :
    
    
    def get_voisin(self, atom, rlim = 3):
        #if type(n_voisin) == list :
        #    n_voisin = np.array(n_voisin)
        if type(atom) == int or type(atom) == np.int64:
            index = atom
        else :
            at= atom.split('_')[0]
            number = int(atom.split('_')[1])
            print(at, number)
            tmp_bol1 = self.data['At'] == at
            tmp = self.data.loc[tmp_bol1]
            tmp_bol2 = tmp['number']  == number
            index = tmp.loc[tmp_bol2].index[0]
        
        x0,y0,z0 = self.data.loc[index][['x','y','z']]
        tmp = self.data.drop(index)
        ax, ay, az = self.diag_lattice
        tmp['x-x0'] = ((tmp['x']-x0)+0.5) % 1-0.5
        tmp['y-y0'] = ((tmp['y']-y0)+0.5) % 1-0.5
        tmp['z-z0'] = ((tmp['z']-z0)+0.5) % 1-0.5
        
        tmp['r(A)'] = np.sqrt((tmp['x-x0']*ax)**2 + (tmp['y-y0']*ay)**2 +(tmp['z-z0']*az)**2)
        """
        res = tmp.sort_values('r(A)')
        index2 = res.index[np.unique(res['r(A)'], return_index=True)[1][n_voisin-1]]
        tmp = np.unique(res['r(A)'])
        print(tmp)
        #r , index, count = tmp[0][n_voisin-1], tmp[1][n_voisin-1], tmp[2][n_voisin-1]
        index = np.array([])
        for r in tmp :
            index = np.append(index, res.loc[res['r(A)'] == r].index)
        index = np.int0(index)
        
        return res.loc[index]
        """
        new_index = []
        for i in tmp.index :
            if tmp['r(A)'][i] < rlim :
                #print(i, tmp['r(A)'][i])
                new_index.append(i)
        return new_index
    
    def get_order(self,path, FILENAME = 'OUTCAR', format_file= 'OUTCAR'):
        if format_file == 'OUTCAR':
            order = []
            for i in grep(LABEL_ion_type,path + FILENAME)[0]:
                order.append(i.split(':')[0].split("=")[-1].split()[0])
        elif format_file == 'POSCAR' or format_file == 'CONTCAR' :
            order = get_row(6,path + FILENAME).split()
        return order
    
    
    
    
    def vasp(self, name = None, decimals = 16):
        if name == None :
            name = ''
            for i in range(len(self.atoms)):
                at = self.atoms[i]
                n = self.n_at[i]
                if n == 1:
                    name += at
                else :
                    name += at+str(n)
        string = name+'\n'
        string += "   %.16f\n" %(1)
        for i in range(3):
            for j in range(3):
                string += " % 4.16f" %(self.lattice[i][j])
            string += '\n'
        for at in self.atoms :
            string += '   %2s' %(at)
        string += '\n'
        for n in self.n_at :
            string += ' %4d' %(n)
        string += '\nDirect\n'
        for i in range(1,len(self)+1) :
            for ax in ['x','y','z'] :
                string += ' % 3.16f' %(np.round(self.data[ax][i], decimals = decimals))
            if 'SD_x' in self.data.columns :
                string += '   %s  %s  %s' %(self.data['SD_x'][i], self.data['SD_y'][i], self.data['SD_z'][i])
            string += '\n'
        return string
        
    
    def write(self,path, FILENAME = 'POSCAR', name = 'POSCAR'):
        #fp = file(path+FILENAME,'w')
        string = name+'\n'
        string += "   %.8f\n" %(1)
        for i in range(3):
            for j in range(3):
                string += " % 4.8f" %(self.lattice[i][j])
            string += '\n'
        for at in self.atoms :
            string += '   %2s' %(at)
        string += '\n'
        string
        #fp.close()
        
    
    def get_n_at_vasp(self,path, FILENAME = 'OUTCAR', format_file= 'OUTCAR'):
        if format_file == 'OUTCAR':
            n_at = np.int0(grep(LABEL_nion,path + FILENAME)[0][0].split("=")[-1].split())
        elif format_file == 'POSCAR' or format_file == 'CONTCAR' :
            n_at = np.int0(get_row(7,path + FILENAME).split())
        return n_at
    
    def get_lattice(self,path, FILENAME = 'OUTCAR', format_file= 'OUTCAR',Cond_pos_at=True,iteration = -1):
        if format_file == 'OUTCAR':
            if Cond_pos_at :
                try : 
                    a_x = np.float64(grep("A1",path + FILENAME)[0][0].split("(")[-1].split(")")[0].split(','))
                    a_y = np.float64(grep("A2",path + FILENAME)[0][0].split("(")[-1].split(")")[0].split(','))
                    a_z = np.float64(grep("A3",path + FILENAME)[0][0].split("(")[-1].split(")")[0].split(','))
                except :
                    n_lattice = grep("direct lattice vectors",path + FILENAME)[1][0]
                    a_x = np.float64(get_row(n_lattice+1,path+FILENAME).split())[0:3]
                    a_y = np.float64(get_row(n_lattice+2,path+FILENAME).split())[0:3]
                    a_z = np.float64(get_row(n_lattice+3,path+FILENAME).split())[0:3]
            else :
                n_lattice = grep("direct lattice vectors",path + FILENAME)[1][0]
                a_x = np.float64(get_row(n_lattice+1,path+FILENAME).split())[0:3]
                a_y = np.float64(get_row(n_lattice+2,path+FILENAME).split())[0:3]
                a_z = np.float64(get_row(n_lattice+3,path+FILENAME).split())[0:3]
        elif format_file == 'POSCAR' or format_file == 'CONTCAR' :
            a_x = np.float64(get_row(3,path + FILENAME).split())
            a_y = np.float64(get_row(4,path + FILENAME).split())
            a_z = np.float64(get_row(5,path + FILENAME).split())
        return np.array((a_x,a_y,a_z)), np.array((a_x[0],a_y[1],a_z[2]))
            
    
    def get_pos_at_outcar(self,path, FILENAME = 'OUTCAR'):
        NSW = int(grep("NSW",path + FILENAME)[0][0].split()[2])
        if NSW == 0 :
            self.iteration = 0
        if self.iteration == 0 :
            l_position = grep(LABEL_position,path + FILENAME)[1][-1]
        else :
            l_position = grep('POSITION',path + FILENAME)[1][self.iteration]+1
        n_tot = self.n_at.sum()
        index = []
        for i in range(0,len(self.order)) :
            index += [self.order[i]]*self.n_at[i]
        n_index = range(1,n_tot+1)
        data = pd.read_csv(path + FILENAME, delim_whitespace = True, skiprows = l_position, usecols=[0,1,2],nrows = n_tot, names=['x','y','z'])
        data.index = range(1,n_tot+1)
        data['At'] = index
        data['n_at'] = n_index
        if self.iteration != 0 :
            a_x, a_y, a_z = self.diag_lattice
            data['x'] /= a_x
            data['y'] /= a_y
            data['z'] /= a_z
        return data  

    def get_pos_at_poscar(self,path, FILENAME = 'POSCAR', force = False):
        n_tot = np.sum(self.n_at)
        index = []
        for i in range(0,len(self.order)) :
            index += [self.order[i]]*self.n_at[i]
        n_index = range(1,n_tot+1)   
        if len(grep('Selective',path+FILENAME)[0]) == 0 :
            
            skiprows = 8
            names = ['x','y','z']
        else :
            skiprows = 9
            names = ['x','y','z', 'SD_x', 'SD_y', 'SD_z']
        if force == True :
            
            names = ['x','y','z', 'SD_x', 'SD_y', 'SD_z']
        data = pd.read_csv(path+FILENAME,names=names,skiprows=skiprows, nrows = n_tot, delim_whitespace = True)     
        data.index = range(1,n_tot+1)
        data['At'] = index
        data['n_at'] = n_index
        return data
        
    def search_at(self, at):
        return self.loc[self['At'] == at]
    #format_file = 'OUTCAR' / 'POSCAR'
    
    def get_atoms_order(self):
        atoms = []
        for at in self.data['At']:
            if at not in atoms :
                atoms.append(at)
        self.atoms = atoms
    
    def get_at_num(self):
        new_index = np.array([])
        tmp = []
        for at in self.atoms :
            if at not in tmp :
                tmp.append(at)
                index = np.int0(self.data.loc[self.data['At'] == at].index)
                new_index = np.append(new_index,index-index[0]+1)
        self.data['num'] = np.int0(new_index)
    
    def get_n_at(self):
        n_at = []
        for at in self.atoms :
            l = len(self[at])
            n_at.append(l)
        self.n_at = n_at
    
    def drop(self,index):
        self.data = self.data.drop(index)
        self.get_atoms_order()
        self.get_at_num()
        self.get_n_at()
        self.data.index = range(1,len(self)+1)
        
    def add(self, at, coordinates, extra = {}):
        series = {'At' : at, 'x' : coordinates[0], 'y' : coordinates[1], 'z' : coordinates[2]}
        
        self.data = self.data.append()
    
    def cart(self) :
        a,b,c = self.lattice
        X = self.data['x']%(1)
        Y = self.data['y']%(1)
        Z = self.data['z']%(1)
        tmp = copy.deepcopy(self)
        data = np.transpose(np.dot(np.transpose(self.lattice), np.array([X,Y,Z])))
        tmp.data[['x','y','z']] = data
        return tmp
    
    def __init__(self,path, 
                 FILENAME = None,
                 format_file = 'OUTCAR',
                 Cond_pos_at = True,
                 iteration = -1, force = False,
                 **kwargs):        
        #FILENAME = kwargs.get('FILENAME',None)
        #format_file = kwargs.get('format_file','OUTCAR')
        #Cond_pos_at = kwargs.get('Cond_pos_at',True)
        #iteration = kwargs.get('iteration',-1)
        
        if FILENAME == None :
            if format_file == 'OUTCAR' :
                FILENAME = 'OUTCAR'
            elif format_file == 'POSCAR' :
                FILENAME = 'POSCAR'
            elif format_file == 'CONTCAR' :
                FILENAME = 'CONTCAR'                
        path = clear_path(path)
        self.iteration = iteration
        self.order = self.get_order(path,FILENAME,format_file)
        self.n_at = self.get_n_at_vasp(path,FILENAME,format_file)
        self.lattice , self.diag_lattice = self.get_lattice(path, FILENAME,format_file,Cond_pos_at,iteration)
        if format_file == 'OUTCAR' :
            self.data = self.get_pos_at_outcar(path,FILENAME)
        elif format_file == 'POSCAR' or format_file == 'CONTCAR' :
            self.data = self.get_pos_at_poscar(path,FILENAME, force = force)
        self.V = np.dot(self.lattice[0], np.cross(self.lattice[1],self.lattice[2]))
        #self.atoms = np.unique(self.data['At'])
        self.get_atoms_order()
        #self.data['number'] = np.int0(np.ones(len(self.data)))
        """
        new_index = np.array([])
        tot_ind = np.array([])
        tmp = []
        for at in self.order :
            if at not in tmp :
                tmp.append(at)
                index = np.int0(self.data.loc[self.data['At'] == at].index)
                
                new_index = np.append(new_index,index-index[0]+1)
        """
        
        #self.data['num'] = np.int0(new_index)
        #print(new_index)
        #self.data['num'] = np.int0(new_index)
        self.get_at_num()
        self.loc = self.data.loc
    
    def verif_data_info(self):
        for at in self.atoms :
            if at not in np.unique(self.data['At']) :
                self.atoms.remove(at)
        self.n_at = []
        for i,at in enumerate(self.atoms) :
            self.n_at.append(len(np.where(self['At'] == at)[0]))
        self.data.index = np.arange(1,len(self.data.index)+1)
        
            
    
    def __repr__(self):
        return "%s" %str(self.data)
    
    def __getitem__(self, index):
        if type(index) == str :
            if index in self.data.columns :
                return self.data[index]
            else :
                return self.data.loc[self.data['At'] == index]
        else :
            return self.data[index]
    
    def __len__(self):
        return np.sum(self.n_at)
    
    
    def get_count_at(self, index_at):
        at = self.data['At'][index_at]
        pos_order = np.where(np.array(self.order) == at)[0][0]
        length = self.n_at[pos_order]
        if length == 1 :
            return at
        else :
            index = np.where(np.array(self.data.loc[self.data['At'] == at].index)==index_at)[0][0]+1
            return at+str(index)
        
    
    # Renvoi un tableau où chaque ligne donne le numéros des atomes constituant un plan atomique.
    # Deux atomes sont considérés dans le même plan si la distance séparant leur coordonnée z est inférieur à dz (en A)           
    def get_plane(self,**kwargs):
        d = kwargs.get('d',0.5)
        direct_lattice = kwargs.get('direct_lattice',True)
        direction = kwargs.get('direction','z')
        
        a_x,a_y,a_z = self.diag_lattice
        A  = {'x':a_x,'y':a_y,'z':a_z}
        POS = self.data.sort_values(direction)
        POS['index'] = np.arange(0,len(POS))
        POS = POS.set_index('index')
        a = 0
        tab = [[]]
        count = 0
        #tab[count].append(POS['n_at'][0])
        a = POS['n_at'][0]
        if direct_lattice :
            d /= A[direction]
        #print POS
        veref = []
        for i in range(0,len(self.data)) :
            if not i in veref :
                #print POS[direction][i]
                j = i+1
                if i != 0:
                    tab.append([])
                    
                tab[count].append(POS['n_at'][i])
                for j in range(0,len(self.data)) :
                    if j not in veref and j != i :
                        if abs((POS[direction][i] - POS[direction][j])) % A[direction] <= d or abs((POS[direction][i] - POS[direction][j] + 1)) % A[direction] <= d :
                            tab[count].append(POS['n_at'][j])
                            veref.append(j)
                        
                """
                while  abs((POS[direction][i] - POS[direction][j%len(self.data)])) % A[direction] <= d and not j in veref:
                    tab[count].append(POS['n_at'][j])
                    veref.append(j%len(self.data))
                    j = j+1
                j = i-1
                while  (POS[direction][i] - POS[direction][j%len(self.data)]) % A[direction] <= d and not j in veref:
                    tab[count].append(POS['n_at'][j%len(self.data)])
                    veref.append(j%len(self.data))
                    j = j-1
                """
                count = count +1
        for i in range(0,len(tab)):
            tab[i].sort()
        for i in range(0,len(tab)):
            tab[i] = np.array(tab[i])
        return np.array(tab)
    
    def get_dist_plane(self, vac = False,**kwargs):
        direction = kwargs.get('direction', 'z')
        planes = kwargs.get('planes', self.get_plane(**kwargs))
        
        def z_moy(plane):
            z = 0.
            zO = 0.
            cont_O = 0
            cont_C = 0
            zC = 0.
            for i in plane:
                if self.data.set_index('n_at')['At'][i] == 'O' :
                    zO += self.data.set_index('n_at')[direction][i]
                    cont_O += 1
                else :
                    zC += self.data.set_index('n_at')[direction][i]
                    cont_C += 1
                    
                z += self.data.set_index('n_at')[direction][i]
            if cont_O == 0 :
                zO = np.nan
            else :
                zO /= cont_O
            if cont_C == 0:
                zC = np.nan
            else :
                zC /= cont_C
            return z/len(plane), zC, zO
        DATA = pd.DataFrame(columns = ['dzp', 'dzm','dz',
                                       'dzCp','dzCm','dzC',
                                       'dzOp','dzOm','dzO'])
        
        for p in range(len(planes)):
            if p == 0 :
                z,zC,zO = z_moy(planes[p])
                zp,zCp,zOp = z_moy(planes[p+1])
                if vac == True:
                    z_, zC_,zO_  = z_moy(planes[-1])
                    series = {'z' : z, 'dzp': zp-z, 'dzm' : z-z_, 'dz' : (zp-z_)/2.,
                              'zC' : zC, 'dzCp': zCp-zC, 'dzCm' : zC-zC_, 'dzC' : (zCp-zC_)/2.,
                              'zO' : zO, 'dzOp': zCp-zC, 'dzOm' : zO-zO_, 'dzO' : (zOp-zO_)/2.,}
                else :
                    series = {'z' : z, 'dzp': zp-z, 'dzm' : np.nan, 'dz' : zp-z,
                              'zC' : zC, 'dzCp': zCp-zC, 'dzCm' : np.nan, 'dzC' : zCp-zC,
                              'zO' : zO, 'dzOp': zOp-zO, 'dzOm' : np.nan, 'dzO' : zOp-zO}
            elif p == len(planes)-1:
                z_ = z
                zC_ = zC
                zO_ = zO
                
                z = zp
                zC = zCp
                zO = zOp
                if vac == True:
                    zp,zCp,zOp = z_moy(planes[0])
                    series = {'z' : z, 'dzp': zp-z, 'dzm' : z-z_, 'dz' : (zp-z_)/2.,
                              'zC' : zC, 'dzCp': zCp-zC, 'dzCm' : zC-zC_, 'dzC' : (zCp-zC_)/2.,
                              'zO' : zO, 'dzOp': zOp-zO, 'dzOm' : zO-zO_, 'dzO' : (zOp-zO_)/2.}
                else :
                    series = {'z' : z, 'dzp': np.nan, 'dzm' : z-z_, 'dz' : z-z_,
                              'zC' : zC, 'dzCp': np.nan, 'dzCm' : zC-zC_, 'dzC' : zC-zC_,
                              'zO' : zO, 'dzOp': np.nan, 'dzOm' : zO-zO_, 'dzO' : zO-zO_}
            else :
                z_ = z
                zC_ = zC
                zO_ = zO
                
                z = zp    
                zC = zCp
                zO = zOp          
                zp,zCp,zOp = z_moy(planes[p+1])
                series = {'z' : z, 'dzp': zp-z, 'dzm' : z-z_, 'dz' : (zp-z_)/2.,
                          'zC' : zC, 'dzCp': zCp-zC, 'dzCm' : zC-zC_, 'dzC' : (zCp-zC_)/2.,
                          'zO' : zO, 'dzOp': zOp-zO, 'dzOm' : zO-zO_, 'dzO' : (zOp-zO_)/2.}
            DATA = DATA.append(series,ignore_index = True)
        return DATA
    
    def get_dist_co(self, planes = None, **kwargs):
        direction = kwargs.get('direction', 'z')
        if planes == None :
            planes = self.get_plane(**kwargs)
        #planes = kwargs.get('planes', self.get_plane(**kwargs))
        def z_moy(plane, direction):
            z = 0.
            if len(plane) > 0:
                for i in plane:
                    z += self.data.set_index('n_at')[direction][i]
                
                return z/len(plane) 
            else :
                return np.nan
        def class_plane(plane):
            C = []
            C_name = []
            O = []
            for i in plane:
                at = self.data.set_index('n_at')['At'][i]
                if at == 'O':
                    O.append(i)
                else :
                    C.append(i)
                    C_name.append(at)
            return [C,C_name], O
        def multi_z(C,O, direction):
            lo = len(O)
            lc = len(C)
            D = np.zeros(lc*lo)
            for i in range(lc):
                zc = self.data.set_index('n_at')[direction][C[i]]
                for j in range(lo):
                    zo = self.data.set_index('n_at')[direction][O[j]]
                    D[i*lo+j] = zo-zc
            return D
        DATA = pd.DataFrame(columns = ['C','O','dx_CO','dy_CO', 'dz_CO', 'dr_CO','dx','dy','dz','dr'])
        
        for p in planes :
            C,O = class_plane(p)
            series = {'C': C, 'O': O}
            for DIR in ['x','y','z']:
                series['d'+DIR+'_CO'] = multi_z(C[0],O,DIR)
                if len(O) == 0 or len(C) == 0:
                    series['d'+DIR] = np.nan
                else :
                    series['d'+DIR] = z_moy(O,DIR)-z_moy(C[0],DIR)
            
            series['dr_CO'] = series['dx_CO']**2+series['dy_CO']**2+series['dz_CO']**2
            series['dr'] = series['dx']**2+series['dy']**2+series['dz']**2
            DATA = DATA.append(series,ignore_index = True)
        return DATA

        


class CHARGE():
    def __init__(self, path, plane = None):
        path = clear_path(path)
        line = grep("total charge", path+"OUTCAR")[1][-1]
        columns = ["n_at"] + get_row(line+2, path+"OUTCAR").split()[3:]
        n_tot = np.int0(grep(LABEL_nion,path + "OUTCAR")[0][0].split("=")[-1].split()).sum()
        #print columns, n_tot
        self.data = pd.read_csv(path + "OUTCAR", delim_whitespace = True, skiprows = line+3,nrows = n_tot, names=columns)
        last = pd.read_csv(path + "OUTCAR", delim_whitespace = True, skiprows = line+3+n_tot+1,nrows = 1, names=columns)
        
        last['n_at'].loc[0] = 0
        self.data = self.data.append(last, ignore_index = True).sort_values('n_at').reset_index().drop(columns = "index")
        if plane == None :
            poscar = POSCAR(path)
            plane = poscar.get_plane()
        self.byplane_tot = self.rearange(plane)
        
    def rearange(self, plane, orb = 'tot'):
        lmax = 0
        for p in plane :
            if len(p) > lmax :
                lmax = len(p)
        tab = np.zeros((lmax,len(plane)))
        for p in range(len(plane)) :
            
            for ion in range(lmax) :
                if ion < len(plane[p]):
                    tab[ion][p] = self.data[orb][plane[p][ion]]
                else :
                    tab[ion][p] = np.NaN
        return tab

path = "/home/gosteau/Documents/TMP_remi/FeSTO/23_09_2019/SrO/perfect/"
path_record = "/home/gosteau/Documents/These/1ere année/FeSTO/Raw_Data/magmom/"
class MAGMOM():
    
    def write_magmom(self, path, out = '.xlsx', **kwargs):
        #self.data[['n_at','At','num']+list(self.orbs)].to_excel(path, **kwargs)
        columns = ['At','num', 'Magmom','Position (ua)', 'Position (A)']
        sub_columns = ['']+list(self.orbs)+['x','y','z']
        label_col = [0,1]+[2]*len(self.orbs)+[3]*3+[4]*3
        label_sub = [0,0]+list(range(1,len(self.orbs)+4))+list(range(len(self.orbs)+1,len(self.orbs)+4))
        labels = [label_col,label_sub]
        midx = pd.MultiIndex(levels=[columns, sub_columns], labels=labels)
        
        matrix_data = np.array(self.data['At'])
        matrix_data = np.column_stack([matrix_data, self.data['num']])
        for orb in self.orbs :
            matrix_data = np.column_stack([matrix_data, self.data[orb]])
        for pos in ['x','y','z'] :
            matrix_data = np.column_stack([matrix_data, np.append([np.nan],self.poscar[pos])])
        ax, ay, az = self.poscar.diag_lattice
        A = {'x':ax,'y':ay,'z':az}
        for pos in ['x','y','z'] :
            matrix_data = np.column_stack([matrix_data, np.append([np.nan],self.poscar[pos]*A[pos])])
            
            
            
        data = pd.DataFrame(matrix_data,columns=midx)
        
        
        def row_colour(row):
            def_col = 'red'
            color = '#a5a5a5'
            color_tab = ['']*2
            #color_tab += ['border-color:'+def_col+','+ 'background-color:'+color]*len(self.orbs)
            color_tab += ['background-color:'+def_col]*len(self.orbs)
            color_tab += ['']*3
            color_tab += ['background-color:'+def_col]*3
            return color_tab
        
        def border_colour(row):
            def_col = 'black'
            color = '#a5a5a5'
            color_tab = ['border-color:'+def_col]*2
            color_tab += ['border-color:'+def_col]*len(self.orbs)
            color_tab += ['border-color:'+def_col]*3
            color_tab += ['border-color:'+def_col]*3
            return color_tab
        
        data = data.style.apply(row_colour,axis=1)
        #data = data.style.apply(border_colour,axis=1)
        #data = data.style.set_properties(**{'background-color': 'black',
        #                    'color': 'grey',
        #                    'border-color': 'white'})
        
        data.to_excel(path+out, **kwargs)
        try :
            #data.to_excel(path+out, **kwargs)
            return 1
        except :
            return 0
        
    def write_magmom(self, path, out = '.xlsx', COLORS = COLORS, **kwargs):
        #self.data[['n_at','At','num']+list(self.orbs)].to_excel(path, **kwargs)
        columns = ['At','num', 'Magmom','Position (ua)', 'Position (A)', 'plane']
        sub_columns = ['']+list(self.orbs)+['x','y','z']
        label_col = [0,1]+[2]*len(self.orbs)+[3]*3+[4]*3+[5]
        label_sub = [0,0]+list(range(1,len(self.orbs)+4))+list(range(len(self.orbs)+1,len(self.orbs)+4))+[0]
        labels = [label_col,label_sub]
        midx = pd.MultiIndex(levels=[columns, sub_columns], labels=labels)
        
        matrix_data = np.array(self.data['At'])
        matrix_data = np.column_stack([matrix_data, self.data['num']])
        for orb in self.orbs :
            matrix_data = np.column_stack([matrix_data, self.data[orb]])
        for pos in ['x','y','z'] :
            matrix_data = np.column_stack([matrix_data, np.append([np.nan],self.poscar[pos])])
        ax, ay, az = self.poscar.diag_lattice
        A = {'x':ax,'y':ay,'z':az}
        for pos in ['x','y','z'] :
            matrix_data = np.column_stack([matrix_data, np.append([np.nan],self.poscar[pos]*A[pos])])
        planes = self.poscar.get_plane(**kwargs)
        index_plane = []
        for i in self.poscar.data.index:
            for p in range(len(planes)) :
                if i in planes[p] :
                    index_plane.append(p)
        index_plane = np.array(index_plane)
        matrix_data = np.column_stack([matrix_data, np.append([np.nan],index_plane)])
            
            
        data = pd.DataFrame(matrix_data,columns=midx)
        
        writer = pd.ExcelWriter(path+out, engine='xlsxwriter')
        data.to_excel(writer, sheet_name='Sheet1')
        
        workbook  = writer.book
        worksheet = writer.sheets['Sheet1']
        
        
        ini = 2+3
        end = len(data)+ini-2
        FORMATS = lambda at : workbook.add_format({'font_color' : COLORS[at]})
        
        FORMATS_BOTT = lambda at : workbook.add_format({'top' : 1,'top_color' : COLORS[at]})
        col = 'B'
        for at in self.poscar.order :
            i1 = self.poscar.loc[self.poscar['At'] == at].index[0]-1+ini
            e1 = self.poscar.loc[self.poscar['At'] == at].index[-1]-1+ini
            cols = "%s%d:%s%d" %(col, i1, col, e1)
            worksheet.conditional_format(cols, {'type':     'no_errors',
                                         'format': FORMATS(at)})
            worksheet.set_row(e1, cell_format = FORMATS_BOTT(at))
            
        
        float_format = workbook.add_format({'num_format': '#0.000'})
        #plane_format = workbook.add_format({'3_color_scale': '#0.000'})
        worksheet.set_column('D:M', None, float_format)
        
        col_magmom_1 = 'D'
        col_magmom_2 = 'G'
        cols = "%s%d:%s%d" %(col_magmom_1, ini, col_magmom_2, end)
        worksheet.conditional_format(cols, {'type':     'data_bar',
                                         'bar_direction': 'center', 
                                         'bar_axis_position':'middle',
                                         'bar_solid':True})
        col_plane = 'N'       
        cols = "%s%d:%s%d" %(col_plane, ini, col_plane, end)
        worksheet.conditional_format(cols, {'type':     '3_color_scale'})
        
        format_def = workbook.add_format({'border': 1,'border_color':'#c0c0c0', 'bg_color' : '#9d9d9d'})
        col_num = 'C'
        cols = "%s%d:%s%d" %(col_num, ini, col_num, end)
        worksheet.conditional_format(cols, {'type':     'no_errors',
                                         'format': format_def})
        col_pos_1 = 'H'
        col_pos_2 = 'J'
        cols = "%s%d:%s%d" %(col_pos_1, ini, col_pos_2, end)
        worksheet.conditional_format(cols, {'type':     'no_errors',
                                         'format': format_def})
                    
                    
                   
        
            
                                                   
        writer.save()
    
    def moy_plane_atom(self):
        res= []
        i = 0
        for p in self.plane :
            Atoms, Count = np.unique(self.poscar.data.loc[self.plane[i]]['At'], return_counts=True)
            res.append({})
            for at in Atoms :
                res[i][at] = {}
                for l in self.orbs :
                    res[i][at][l] = 0
            for n_at in p:
                at = self.poscar.data['At'][n_at]
                for l in self.orbs :
                    #print(i,at,l)
                    res[i][at][l] += self.data[l][n_at]
            for j in range(len(Atoms)):
                at = Atoms[j]
                c = Count[j]
                for l in self.orbs :
                    res[i][at][l] /= c
            i += 1
        return res
    
    def moy_per_at(self):
        res = {}
        for at in self.poscar.order :
            res[at] = {}
            for l in self.orbs:
                res[at][l] = np.full(len(self.plane), np.nan)
                for i in range(len(self.plane)):
                    if at in self.moy[i].keys():
                        res[at][l][i] = self.moy[i][at][l]
        return res    
    
    
    def get_atomic_magmom_from_plane(self, used_orbs = None):
        if used_orbs == None :
            used_orbs = self.orbs
        MAX = 0
        for i in self.plane:
            if len(i)>MAX :
                MAX = len(i)
        L = len(used_orbs)
        TAB = np.full((len(self.plane), MAX*(L+1)), np.nan, dtype='<U21')
        
        for i in range(len(self.plane)):
            k = 0
            for j in self.plane[i]:
                at = self.poscar['At'][j]+'_'+str(self.poscar['num'][j])
                number = self.poscar['num'][j]
                #print(i,j,k,at,number,end=' ')
                TAB[i][(L+1)*k]=at
                for l in range(L) :
                    orb = used_orbs[l]
                    #print(orb)
                    mom = self.data[orb][j]
                    TAB[i][(L+1)*k+1+l]=mom
                    #print(l,mom,end = ' ')
                #print('')
                k+= 1
        return TAB, MAX
    
    def get_multiindex_tab_for_details(self, used_orbs = None):
        if used_orbs == None :
            used_orbs = self.orbs
            
            
        DATA = np.arange(len(self.plane))
        col = ['plane']
        sub_col = used_orbs
        sub_col = np.append(sub_col,'name')
        labels = [[0],[-1]]
        i = 1
        L = len(used_orbs)
        L_orb = len(self.orbs)
        for at in self.poscar.order :
            col += ['moy_'+at]
            labels[0] += [i]*(L)# #[i,i,i,i]
            labels[1] += range(L)#[0]#[0,1,2,3]
            for orb in sub_col[0:L] :
                DATA = np.column_stack((DATA, self.res[at][orb]))
            i+=1
        
        #sub_col += ['name']
        for j in range(self.MAX):
            col += ['Atom_'+str(j+1)]
            labels[0] += [i]*(L+1)#[i,i]#[i,i,i,i,i]
            labels[1] += [L]+list(range(L))#[1,0]#[4,0,1,2,3]
            i+=1
        
        m = [0]
        for orb in used_orbs :
            #print(np.where(self.orbs == orb)[0][0])
            m.append(np.where(self.orbs == orb)[0][0]+1)
        m = np.array(m)
        orbs = np.array([])
        for i in range(self.MAX):
            #print(np.array([1]*(L+1)), np.array([L_orb+1]*(L+1))*i, i*np.array([L_orb+1]*(L+1))+m)
            orbs = np.append(orbs, i*np.array([L_orb+1]*(L+1))+m)
        #print(orbs)
        DATA = np.column_stack((DATA, self.TAB[:,np.int0(orbs)]))
        #print(used_orbs, self.orbs, m, np.shape(DATA))
        
        #print(col, sub_col)
        #print(labels[0])
        #print(labels[1])
        #print(self.used_orbs)
        midx = pd.MultiIndex(levels=[col, sub_col], labels=labels)
        df = pd.DataFrame(DATA,columns=midx)
        return df
    
    def get_data(self):
        
        columns = ['At','num', 'Magmom','Position (ua)', 'Position (A)', 'plane']
        sub_columns = ['']+list(self.orbs)+['x','y','z']
        sub_sub_columns = ['','x','y','z']
        label_col = [0,1]+[2]*len(self.orbs)*3+[3]*3+[4]*3+[5]
        sub = []
        for i in range(len(list(self.orbs))) :
            sub.append([i+1]*3)
        label_sub = [0,0]+sub+[0]
        labels = [label_col,label_sub]
        midx = pd.MultiIndex(levels=[columns, sub_columns], labels=labels)
    
    
    def __init__(self, path, **kwargs):
        path = clear_path(path)
        self.poscar = kwargs.get('poscar',POSCAR(path, **kwargs))
        self.plane = kwargs.get('plane',self.poscar.get_plane(**kwargs))
        
        line = grep("magnetization (x)", path+"OUTCAR")[1][-1]
        columns = ["n_at"] + get_row(line+2, path+"OUTCAR").split()[3:]
        n_tot = np.int0(grep(LABEL_nion,path + "OUTCAR")[0][0].split("=")[-1].split()).sum()
        self.data = pd.read_csv(path + "OUTCAR", delim_whitespace = True, skiprows = line+3,nrows = n_tot, names=columns)
        last = pd.read_csv(path + "OUTCAR", delim_whitespace = True, skiprows = line+3+n_tot+1,nrows = 1, names=columns)
        self.orbs = self.data.columns[1:]
        last['n_at'].loc[0] = 0
        self.data = self.data.append(last, ignore_index = True).sort_values('n_at').reset_index().drop(columns = "index")
        self.data['At'] = self.poscar['At']
        self.data['num'] = np.int0([0]+list(self.poscar['num']))
        self.moy = self.moy_plane_atom()
        self.res = self.moy_per_at()
        self.TAB, self.MAX = self.get_atomic_magmom_from_plane()
        

class MAGMOM2():
    
    def get_lines(self, skiprows):
        count = 0
        with open(self.path+self.file, 'r') as f:
            txt = ''
            while count < skiprows :
                line = f.readline()
                count += 1
            for n in range(np.sum(self.poscar.n_at)):
                line = f.readline()
                txt += line
            line = f.readline()
            line = f.readline()
            txt += line
        return txt
    
    def get_mag(self, magnetization = "magnetization"):
        n_line = grep(magnetization, self.path+self.file)[1][-1]
        self.orbs = np.array(get_row(n_line+2, self.path+self.file).split())[3:]
        natoms = np.sum(self.poscar.n_at)
        txt = self.get_lines(n_line+3)
        tmp = np.array(txt.split()).reshape((natoms+1,1+len(self.orbs)))
        data = np.float64(tmp[:,1:])
        return data             
    
    def get_data_cart(self):
        data = np.transpose(np.dot(self.SC,np.transpose(self.data[:,0:0+3])))
        tot = np.dot(self.SC,np.transpose(self.tot[:3]))
        
        for i in range(1,len(self.orbs)) :
            tmp = np.transpose(np.dot(self.SC, np.transpose(self.data[:,i*3:i*3+3]) ))
            data = np.column_stack((data, tmp))
            tmp = np.dot(self.SC,np.transpose(self.tot[i*3:i*3+3]) )
            tot = np.append(tot, tmp )                        
        return data, tot       
    
    def get_data_lattice(self):
        datac, totc = self.get_data_cart()
        lattice = self.poscar.lattice
        lattice[0] /= np.linalg.norm(lattice[0])
        lattice[1] /= np.linalg.norm(lattice[1])
        lattice[2] /= np.linalg.norm(lattice[2])
        
        Mlat = np.linalg.inv(np.transpose(lattice))
        V = self.poscar.V
        V = 1
        
        data = np.transpose(np.dot(Mlat,np.transpose(datac[:,0:0+3])))
        tot = np.dot(Mlat,np.transpose(totc[:3]))
        
        for i in range(1,len(self.orbs)) :
            tmp = np.transpose(np.dot(Mlat, np.transpose(datac[:,i*3:i*3+3]) ))
            data = np.column_stack((data, tmp))
            tmp = np.dot(Mlat,np.transpose(totc[i*3:i*3+3]) )
            tot = np.append(tot, tmp )                        
        return data, tot      
        
        
                                  
    
    def get_data(self, data = None, tot = None):
        if type(data) == type(None) :
            data = self.data
        if type(tot) == type(None) :
            tot = self.tot
        columns = ['At','num', 'Magmom','Position (ua)', 'Position (A)']
        sub_columns = ['']+list(self.orbs)+['x','y','z']
        sub_sub_columns = ['','x','y','z']
        label_col = [0,1]+[2]*len(self.orbs)*3+[3]*3+[4]*3
        sub = []
        for i in range(len(list(self.orbs))) :
            sub += [i+1]*3
        label_sub = [0,0]+sub+[0]*(3+3)
        label_sub_sub = [0,0] + [1,2,3] * len(self.orbs) + [0]*(3+3)
        
        labels = [label_col,label_sub, label_sub_sub]
        midx = pd.MultiIndex(levels=[columns, sub_columns, sub_sub_columns], labels=labels)
        px = np.array(self.poscar['x'])* np.linalg.norm(self.poscar.lattice[0])
        py = np.array(self.poscar['y'])* np.linalg.norm(self.poscar.lattice[0])
        pz = np.array(self.poscar['z'])* np.linalg.norm(self.poscar.lattice[0])
        new = np.transpose(np.array([px,py,pz]))
        TOT = np.array([['',0]+ list(tot)+ [0]*6])
        DATA = np.column_stack((self.poscar['At'], self.poscar['num'], data, self.poscar[['x','y','z']],new))
        FULL = np.append(TOT, DATA, axis = 0)
        
        tab = pd.DataFrame(FULL, columns=midx)
        return tab
              
    def __init__(self, path, file = 'OUTCAR', **kwargs):
        path = clear_path(path)
        self.file = file
        self.path = path
        self.poscar = kwargs.get('poscar',POSCAR(path, **kwargs))
        
        if get_lsorbit(path+'OUTCAR') :
            tmp = self.get_mag("magnetization (x)")
            self.x = tmp[:-1]
            self.totx = np.append(tmp[-1][0],np.float64(tmp[-1][1:]))
            tmp = self.get_mag("magnetization (y)")
            self.y = tmp[:-1]
            self.toty = np.append(tmp[-1][0],np.float64(tmp[-1][1:]))
            tmp = self.get_mag("magnetization (z)")
            self.z = tmp[:-1]
            self.totz = np.append(tmp[-1][0],np.float64(tmp[-1][1:]))
           
            
            #data = np.array(list(range(len(self.x))))
            #tot = np.array(['tot'])
            tmpx = self.x[:,0]
            tmpy = self.y[:,0]
            tmpz = self.z[:,0]
            data = np.transpose([tmpx,tmpy,tmpz])
            
            tmpx = self.totx[0]
            tmpy = self.toty[0]
            tmpz = self.totz[0]
            tot = np.transpose([tmpx,tmpy,tmpz])
            for orb in range(1,len(self.orbs)):
                tmpx = self.x[:,orb]
                tmpy = self.y[:,orb]
                tmpz = self.z[:,orb]
                data = np.column_stack((data, tmpx,tmpy,tmpz))
                
                tmpx = self.totx[orb]
                tmpy = self.toty[orb]
                tmpz = self.totz[orb]
                tot = np.append(tot, [tmpx,tmpy,tmpz])
            self.data = data
            self.tot = tot
            self.tab = self.get_data()
            self.SC, self.CS = get_maxis(self.path+self.file)
            
            #self.data = self.data.append(last, ignore_index = True).sort_values('n_at').reset_index().drop(columns = "index")
            #self.data['At'] = self.poscar['At']
            #self.data['num'] = np.int0([0]+list(self.poscar['num']))
            
    def __repr__(self):
        return "%s" %str(self.tab)
"""
class POLA():
    def get_pola(self):
        outcar = open(self.path+'OUTCAR','r')
        for line in outcar :
            if "p[ion]" in line :
                self.pion = np.float64(line.split('(')[1].split(')')[0].split())
            elif "p[elc]" in line :
                self.pelc = np.float64(line.split('(')[1].split(')')[0].split())
        outcar.close()
        self.ptot = self.pion + self.pelc
        
        
    
    
    def __init__(self, path):
        self.path = path
        self.get_pola()
        poscar = POSCAR(path)
        a,b,c = poscar.lattice.transpose()
        self.V = np.dot(a,np.cross(b,c))
        out = INFO_OUTCAR(path)
        self.modulos = np.array([np.linalg.norm(a), np.linalg.norm(b), np.linalg.norm(c)])
        #self.modulos = 2*np.pi/out.length_rec_vectors
        self.pion_mod = self.pion % self.modulos
        self.pelc_mod = self.pelc % self.modulos
        self.ptot_mod = self.ptot % self.modulos
"""
        

        
def onpick_pt(event):
    thisline = event.artist
    ldata = thisline.get_label()
    sp = ldata.split()
    label = sp[0]
    x,y,z = np.float64(sp[1:])
    print("\n%8s %6s %6s %6s\n%8s %.4f %.4f %.4f" %('label','x','y', 'z',label,x,y,z))
#fig.canvas.mpl_connect('pick_event', onpick) 
#def onpick_pt(event):
#    print(event.artist.get_xdata())

def display_pos_at_plane(poscar, p, magmom = None, ax = None,auto = True,axis = 'xy',**kwargs):
    AT_displays = {'color'       : {'Sr' : 'green', 'Ti' : 'blue', 'Fe' : 'orange', 'O' : 'red'},
                   'marker'      : {'Sr' : 'o'    , 'Ti' : 'o'   , 'Fe' : 'o'     , 'O' : 'o'},
                   'marker_size' : {'Sr' : 20      , 'Ti' : 20     , 'Fe' : 20       , 'O' : 10}
    }
    planes = poscar.get_plane(**kwargs)
    poscar_plane = poscar.loc[planes[p]]
    tick_size = kwargs.get('tick_size',12)
    tick_width = kwargs.get('tick_width',1.5)
    if ax == None : 
        ax = plt.gca()
        
        ax.set_aspect('equal', adjustable='box')
        ax.minorticks_on()
        ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = tick_size, width = tick_width)
        ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = tick_size/2)
        ax.locator_params(axis='both', nbins=5)
    if auto == True:
        xmax = np.max(poscar['x'])
        ymax = np.max(poscar['y'])
        zmax = np.max(poscar_plane['z'])
        
        xmin = np.min(poscar['x'])
        ymin = np.min(poscar['y'])
        zmin = np.min(poscar_plane['z'])
        dx = (xmax-xmin) / 10
        dy = (ymax-ymin) / 10
        dz = (zmax-zmin) / 10
        if axis == 'xy' :
            ax.set_xlim(xmin-dx,xmax+dx)
            ax.set_ylim(ymin-dy,ymax+dy)
        elif axis == 'xz' :
            ax.set_xlim(xmin-dx,xmax+dx)
            ax.set_ylim(zmin-dy,zmax+dy)
        elif axis == 'yz' :
            ax.set_xlim(ymin-dx,ymax+dx)
            ax.set_ylim(zmin-dy,zmax+dy)
            
        ax.grid(True)
    
    for i in poscar_plane.index :
        at = poscar_plane['At'][i]
        index = np.where(poscar.loc[poscar['At'] == at].index == i)[0][0]
        color = AT_displays['color'][at]
        marker = AT_displays['marker'][at]
        marker_size = AT_displays['marker_size'][at]
        x = poscar_plane['x'][i]
        y = poscar_plane['y'][i]
        z = poscar_plane['z'][i]
        label = "%s_%d %.4f %.4f %.4f"%(at,index,x,y,z)
        if axis == 'xy' :
            ax.plot(x,y,markersize = marker_size, color = color, marker = marker,markeredgecolor = 'black', picker = 10, label= label)

        if axis == 'xz' :
            ax.plot(x,z,markersize = marker_size, color = color, marker = marker,markeredgecolor = 'black', picker = 10, label= label)

        if axis == 'yz' :
            ax.plot(y,z,markersize = marker_size, color = color, marker = marker,markeredgecolor = 'black', picker = 10, label= label)

    
    
def get_IF(n_at,PLANE):
    c = 0
    if type(n_at) == int or type(n_at) == np.int64:
        for p in PLANE:
            if n_at in p :
                return c
            else :
                c+= 1
    else :
        c = []
        for p in n_at :
            c.append(get_IF(p,PLANE))
        return c
    return None


def display_vline(ax, poscar, atom, pos, plane = None,**kwargs):
    ylim = ax.get_ylim()
    if plane == None :
        plane = poscar.get_plane()
    if type(pos) == int :
        pos = [pos]
    n_ats = np.array(poscar.data['n_at'][poscar.data['At'] == atom])[pos]
    x = get_IF(n_ats, plane)
    plt.vlines(x, *ylim, **kwargs)
    
    

def display_inter_dist(path = None, poscar = None, dist_plane = None, poscar_ref = None,dist_plane_ref = None, plane_range = None,ax = None, what = 'dz', **kwargs):
    decal = kwargs.get('decal',0)
    kwargs.pop('decal',None)
    if type(dist_plane) == type(None) :
        if poscar == None :
            poscar = POSCAR(path)
        dist_plane = poscar.get_dist_plane()
    if plane_range == None :
        plane_range = dist_plane.index
    if ax == None :
        ax = plt.gca()
    if len(what.split('+'))  > 1 :
        what = what.split('+')
        op = 1
    elif len(what.split('-'))  > 1 :
        what = what.split('-')
        op = -1
    else :
        what = [what,what]
        op = 0
    if type(dist_plane_ref) == type(None) :
        Y = dist_plane[what[0]].loc[plane_range]*poscar.diag_lattice[-1]+op*dist_plane[what[1]].loc[plane_range]*poscar.diag_lattice[-1]
        X = np.array(plane_range)+decal
        ax.plot(X, Y, **kwargs)
    else :
        Y1 = dist_plane[what[0]].loc[plane_range]*poscar.diag_lattice[-1]+op*dist_plane[what[1]].loc[plane_range]*poscar.diag_lattice[-1]
        Yref = dist_plane_ref[what[0]].loc[plane_range]*poscar.diag_lattice[-1]+op*dist_plane_ref[what[1]].loc[plane_range]*poscar.diag_lattice[-1]
        Y = Y1-Yref
        X = np.array(plane_range)+decal
        ax.plot(X, Y, **kwargs)
    #plt.xlim(*plt.xlim())
    #plt.ylim(*plt.ylim())
    return ax

def routine_display_inter_dist_per_ref(path, path_ref = None, planes = None, ax = None, IF = None,**kwargs):
    tick_size = kwargs.get('tick_size',12)
    tick_width = kwargs.get('tick_width',1.5)
    kwargs.pop('tick_size',None)
    kwargs.pop('tick_width',None)
    poscar = POSCAR(path)
    dist_plane = poscar.get_dist_plane(**kwargs)
    if path_ref == None :
        dist_plane_ref = None
        poscar_ref = None
    else :
        poscar_ref = POSCAR(path_ref)
        dist_plane_ref = poscar_ref.get_dist_plane(**kwargs)
    if planes == None :
        planes = list(dist_plane.index)
    if ax == None :
        ax = plt.gca()
    ax.minorticks_on()
    ax.tick_params(which = 'major',top=True, right=True, direction = 'inout', size = tick_size, width = tick_width)
    ax.tick_params(which = 'minor',top=True, right=True, direction = 'in', size = tick_size/2)
    ax.hlines(0,planes[0], planes[-1], linestyle = ':')
    ax.locator_params(axis='x', nbins=1)
    kwargs.pop('d',None)
    what = kwargs.get('what', 'dz')
    
    
    kwargs.pop('what',None)
    display_inter_dist(poscar = poscar, poscar_ref = poscar_ref, 
                       dist_plane=dist_plane, dist_plane_ref=dist_plane_ref,
                       plane_range=planes, what = what, **kwargs)
    ax.set_xlim(planes[0], planes[-1])
    #ax.set_ylim(ax.get_ylim())
    if IF != None :
        colorsbox = {'box' : [[planes[0],IF],[IF,planes[-1]]], 'color' : ['green','lightgrey'], 'alpha' : 0.2}
        for i in range(0,len(colorsbox['box'])):
            x0 = colorsbox['box'][i][0]
            dx = colorsbox['box'][i][1]-x0
            y0 = ax.get_ylim()[0]
            dy = ax.get_ylim()[1]-ax.get_ylim()[0]
            color = colorsbox['color'][i]
            alpha = colorsbox['alpha']
            rec =plt.Rectangle((x0,y0),dx,dy,color=color,alpha=alpha)
            ax.add_patch(rec)
        tab = np.sort(np.array(planes)[::2]+decal)
        label = (IF-np.sort(np.array(planes)[::2]))/2
        #plt.xlim(0, len(plane)-1)
        #for i in tab:
        #for i in range_plane :
        #    label += [i/2]
        
        plt.xticks(tab, label)




def display_CO(path = None, poscar = None, CO = None, plane_range = None,ax = None, sym = False, what = 'dz', **kwargs):
    if type(CO) == type(None) :
        if poscar == None :
            poscar = POSCAR(path)
        CO = poscar.get_dist_co()
    if plane_range == None :
        plane_range = CO.index
    if ax == None :
        ax = plt.subplot()
    if sym :
        plane_range = CO.index[:len(CO)/2+1]
        plt.plot(plane_range, CO[what][plane_range]*poscar.diag_lattice[-1], **kwargs)
        plane_range = CO.index[len(CO)/2:]
        plt.plot(plane_range, -CO[what][plane_range]*poscar.diag_lattice[-1], **kwargs)
    else :
        plt.plot(plane_range, CO[what][plane_range]*poscar.diag_lattice[-1], **kwargs)
    plt.xlim(*plt.xlim())
    plt.ylim(*plt.ylim())
    return ax

def routine(ax, IF1, IF2, plane, di = 0, xlim = None, ylim = None):
    majorLocator = IndexLocator(base = 2, offset = di) #MultipleLocator(2)
    minorLocator = MultipleLocator(1)
    majorFormatter = FormatStrFormatter('%d')
    if xlim == None :
        plt.xlim(0,len(plane))
    else : 
        plt.xlim(*xlim)
    
    #majorLocator = MultipleLocator(2)
    #minorLocator = MultipleLocator(1)
    #majorFormatter = FormatStrFormatter('%d')
    #ax.set_xticks(range(len(plane)))
    ax.xaxis.set_major_locator(majorLocator)
    ax.xaxis.set_major_formatter(majorFormatter)
    ax.xaxis.set_minor_locator(minorLocator)
    if ylim == None:
        ylim = plt.ylim()
    tab = range(di,len(plane),2)
    label = []
    #plt.xlim(0, len(plane)-1)
    for i in tab:
    #for i in range_plane :
        if i <= (IF1+IF2)/2 :
            label += [(IF1-i)/2]
        else : 
            label += [(i-IF2)/2]
    if IF2 < len(plane):
        colorsbox = {'box' : [[0,IF1],[IF1,IF2],[IF2,len(plane)-1]], 'color' : ['purple','green','purple'], 'alpha' : 0.2}
    else :
        colorsbox = {'box' : [[0,IF1],[IF1,IF2]], 'color' : ['lightgrey','green'], 'alpha' : 0.2}
    for i in range(0,len(colorsbox['box'])):
        x0 = colorsbox['box'][i][0]
        dx = colorsbox['box'][i][1]-x0
        y0 = ylim[0]
        dy = ylim[1]-ylim[0]
        color = colorsbox['color'][i]
        alpha = colorsbox['alpha']
        rec =plt.Rectangle((x0,y0),dx,dy,color=color,alpha=alpha)
        plt.gca().add_patch(rec)
    plt.xticks(tab,label)
    ax2 = ax.twiny()
    ax2.get_xaxis().set_visible(False)
    
    ax2.set_xlim(*ax.get_xlim())

def display_magmom_moy(path, ax = None,mag = None,marker= 'o', plane_range = None, **kwargs):
    if ax == None :
        ax = plt.subplot()
    if mag == None :
        mag =  MAGMOM(path, **kwargs)
    if plane_range == None :
        plane_range = range(len(mag.plane))
    line = kwargs.get('line', True)
    l = kwargs.get('l','tot')
    
    for at in mag.poscar.order :
        color = COLORS[at]
        if line :
            ax.plot(range(len(mag.plane)), mag.res[at][l],color = color, marker = marker,**kwargs)
        else :
            ax.scatter(range(len(mag.plane)), mag.res[at][l],color = color, marker = marker,**kwargs)
    return ax
    
    
        


