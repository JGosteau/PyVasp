from lib_gen import *

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
            self.path = clear_path(path)
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
            self.path = None
            self.nkpt = None
            self.reciprocal = reciprocal
            self.data = pd.DataFrame(columns=['kpt','kx','ky','kz','w_k'])    
        self.loc = self.data.loc
        self.u, self.v , self.w = self.get_uvw(path)

    def get_KPOINTS_file(self):
        with open(self.path+'KPOINTS', encoding='utf-8') as f:
            string = f.readlines()
        INI = np.array([i.split() for i in string[4::3]])
        INI = {'coord' : np.float64(INI[:,[0,1,2]]), 'name' : INI[:,-1]}
        END = np.array([i.split() for i in string[5::3]])
        END = {'coord' : np.float64(END[:,[0,1,2]]), 'name' : END[:,-1]}
        return INI, END
        
    def auto_kpath(self, LINES = None):
        INI, END = self.get_KPOINTS_file()
        n_line = len(INI['coord'])
        nkpts = len(self)
        nkpt_per_line = int(nkpts/n_line)
        if type(LINES) == type(None) :
            LINES = range(n_line)
        KPATHS = [[i*nkpt_per_line,(i+1)*nkpt_per_line-1] for i in range(n_line)]
        kpts, knorm = self.get_knorm(KPATHS)
        POS_INI = np.array(KPATHS)[:,0]
        POS_END = np.array(KPATHS)[:,1]
        POS = np.column_stack((POS_INI, POS_END))
        names = np.column_stack((INI['name'], END['name']))[LINES]
        COORDS = np.column_stack((INI['coord'], END['coord']))[LINES]
        return kpts, knorm, POS, names, COORDS
        
    
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
    