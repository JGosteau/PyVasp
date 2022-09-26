from lib_gen import *
from lib_class_kpt3 import *


class BANDS():
    def __init__(self, path, filename = 'PROCAR', filename_outcar = 'OUTCAR', prec_band = 'OUTCAR'):
        if type(path) == str :
            self.path = clear_path(path)
        else : 
            self.path = [clear_path(p) for p in path]
        self.filename = filename
        info = self.get_info_procar()
        self.nkpts = info['nkpts']
        self.nions = info['nions']
        self.nbands = info['nbands']
        self.orbs = info['orbs']
        self.lsorbit = get_lsorbit(self.path+filename_outcar)
        self.ispin = get_ISPIN(path)
        #self.energy = get_Energy_from_OUTCAR(self.path).reshape((self.nkpts, self.nbands, 3))
        self.energy = self.get_band_energy(prec_band)
        self.kpoints = KPOINTS(self.path,  fileformat = 'OUTCAR')
        self.kpoints_dir = KPOINTS(self.path,  fileformat = 'OUTCAR', reciprocal = False)
    
    def get_band_energy(self, prec_band) :
        if prec_band == 'OUTCAR' :
            if self.ispin == 1 :
                E = get_Energy_from_OUTCAR(self.path).reshape((1,self.nkpts, self.nbands, 3))
            else :
                E = get_Energy_from_OUTCAR(self.path).reshape((2,self.nkpts, self.nbands, 3))
        elif prec_band == 'PROCAR' :
            command = 'grep "^band" %s' %(self.path+self.filename)
            f = os.popen(command)
            sent = f.read()
            f.close()
            if self.ispin == 1 :
                E = np.float64(np.array(sent.split()).reshape(1,self.nkpts, self.nbands, 8)[:,:,:,1::3])
            else :
                E = np.float64(np.array(sent.split()).reshape(2,self.nkpts, self.nbands, 8)[:,:,:,1::3])
        return E
    
    def get_bands_limit(self, Elim, Ef = 0,**kwargs):
        Emin, Emax = Elim
        if self.ispin == 1 :
            ixs1 = np.where(self.energy[0,:,:,1]-Ef <= Emax)
            ixs1 = np.unique(ixs1[1])
            ixs2 = np.where(self.energy[0,:,:,1]-Ef >= Emin)
            ixs2 = np.unique(ixs2[1])
            return np.array(np.intersect1d(ixs1, ixs2))+1
        else :
            ixs1 = np.where(self.energy[0,:,:,1]-Ef <= Emax)
            ixs1 = np.unique(ixs1[1])
            ixs2 = np.where(self.energy[1,:,:,1]-Ef >= Emin)
            ixs2 = np.unique(ixs2[1])
            ix_spin_up = np.array(np.intersect1d(ixs1, ixs2))+1
            ix_spin_dn = np.array(np.intersect1d(ixs1, ixs2))+1
            return np.union1d(ix_spin_up, ix_spin_dn)
            
    
    
    def get_Ef_from_occ(self, w_trig = 0):
        ixs = np.where(self.energy[:,:,:,2] > w_trig)
        return np.max(self.energy[ixs[0],ixs[1],ixs[2],1])
    
    def get_info_procar(self):
        command = "head -2 %s | tail -1" %(self.path+self.filename)
        f = os.popen(command)
        sent = f.read()
        f.close()
        data = np.int64(np.array(sent.split())[3::4])
        
        command = "head -8 %s | tail -1" %(self.path+self.filename)
        f = os.popen(command)
        sent = f.read()
        f.close()
        sent = np.array(sent.split())[0:]
        info = {'nkpts' : data[0], 'nbands' : data[1], 'nions' : data[2], 'orbs' : sent}
        return info
    
    def get_bands(self,BANDS,
                 ATOMS = []):
        if len(BANDS) == 1 :
            keys_bands = BANDS[0]
        else :
            keys_bands = '(%d' %(BANDS[0])
            for b in BANDS[1:] :
                keys_bands += '|%d' %(b)
            keys_bands += ')'
        if self.nions == 1 :
            ATOMS = [1]
        if len(ATOMS) == 0 :
            ATOMS = [str(i) for i in range(1,self.nions+1)] + ['tot']
            keys_atoms = '(%s' %(str(ATOMS[0]))
            for at in ATOMS[1:] :
                keys_atoms += '|%s' %(str(at))
            keys_atoms += ')'
        elif len(ATOMS) == 1 :
            keys_atoms = '%s' %(str(ATOMS[0]))
        else :
            keys_atoms = '(%s' %(str(ATOMS[0]))
            for at in ATOMS[1:] :
                keys_atoms += '|%s' %(str(at))
            keys_atoms += ')'
            
        if self.ispin == 1 and self.lsorbit == True :
            nlines = 2 + (self.nions+1)*4
            SPINS = ['x','y','z','tot']
        else :
            nlines = 2 + (self.nions+1)
            if self.ispin == 1 :
                SPINS = ['tot']
            else :
                SPINS = ['up', 'dn']
        
        command = 'egrep "^band *%s " %s -A %d | egrep "^ *%s "' %(keys_bands, self.path+self.filename, nlines, keys_atoms)
        f = os.popen(command)
        data_string = f.read().replace('tot','0')
        f.close()
        #data_string = np.array(data_string.split()).reshape((len(self.orbs), self.nkpts, len(BANDS), len(ATOMS)) )
        if self.ispin == 1 and self.lsorbit == True :
            data = np.float64(data_string.split()).reshape((self.nkpts, len(BANDS), 4,len(ATOMS), len(self.orbs) ))
        elif self.ispin == 1 : 
            data = np.float64(data_string.split()).reshape((self.nkpts, len(BANDS), 1,len(ATOMS), len(self.orbs) ))
        elif self.ispin == 2 : 
            data = np.float64(data_string.split()).reshape((2,self.nkpts, len(BANDS),len(ATOMS), len(self.orbs) )).transpose(1,2,0,3,4)
        CONTRIBS = {}
        for s in SPINS :
            CONTRIBS[s] = []
        for i, b in enumerate(BANDS) :
            for s in SPINS :
                contrib = CONTRIB(b,s)
                if s == 'dn' :
                    s_i = 1
                else :
                    s_i = 0
                contrib.energy = self.energy[s_i,:,b-1,1]
                contrib.occ = self.energy[s_i,:,b-1,2]
                
                contrib.contrib = {}
                contrib.n_ions = self.nions
                contrib.n_kpts = self.nkpts
                contrib.orbs = self.orbs[1:]
                
                for at in ATOMS :
                    if at == 'tot' or at == 0:
                        tmp = data[:,i, contrib.spin, -1, 1:]
                        contrib.contrib[0] = np.transpose(tmp)
                        contrib.ions.append(0)
                    else :
                        at = int(at)
                        tmp = data[:,i, contrib.spin, at-1, 1:]
                        contrib.contrib[at] = np.transpose(tmp)
                        contrib.ions.append(at)
                
                CONTRIBS[s].append(contrib)
                
        return CONTRIBS
                     
    
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
            for l in liste.replace(',','+').split('+') :
                new.append(self.get_orb(l))
        else :
            for l in liste :
                new.append(self.get_orb(l))
        return new
    
    def test_list_ions(self, liste):
        new = []
        if type(liste) == str :
            if liste == 'tot' :
                new.append(0)
            else :
                for l in liste.replace(',','+').split('+') :
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
                
    
    def __init__(self, band, spin  =""):
        self.band = band
        def string_spin(spin):
            if 'x' in spin :
                return 1
            elif 'y' in spin :
                return 2
            elif 'z' in spin :
                return 3
            elif 'up' in spin :
                return 0
            elif 'dn' in spin :
                return 1
            else :
                return 0
        self.spin = string_spin(spin)
        self.energy = None
        self.occ = None
        self.n_kpts = None
        self.orbs = None
        self.n_ions = None
        self.contrib = None
        self.ions = []
        self.second_access = False
        
        
    
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
        #for i in range(self.n_ions+1):
        for i in self.ions:
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
        #for i in range(1,self.n_ions+1):
        for i in self.ions:
            if type(poscar) == type(None) :
                if i == 0 :
                    string += " tot"
                else :
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
        """
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
                string += "  %1.3f" %(abs(tab[0][j]))"""
            
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
    """    
    def __getitem__(self, key):
        ions = np.array(key[0])
        orbs = np.array(key[1])
        ions = self.test_list_ions(ions)
        orbs = self.test_list_orbs(orbs)
        if len(ions) == 1 and len(orbs) == 1 :
            return self.contrib[ions[0]][orbs[0]]
        else :
            result = np.zeros(self.n_kpts)
            for i in ions:
                for l in orbs :
                    result += self.contrib[i][l]
        return result"""
      
    def __repr__(self) :
        string = ''
        for kpt in range(self.n_kpts):
            string += self.print_tokpt(kpt, return_string = True, return_tab = False)
        return string        
   