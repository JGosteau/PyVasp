from lib_gen import *
from lib_class_kpt3 import *


class BANDS():
    """
    This class permit to obtain the main information of many PROCAR files. It includes also some functions to obtain the projection on orbitals and/or atoms for different bands. The list of the functions are :
    1) addition : "B1 + B2" will return an object BANDS which contains the information of both B1 and B2 if their are compatible (same number of atoms, orbitals, same ISPIN and LSORBIT)
    2) BANDS.get_bands_limit(Elim, Ef) will return an array containing the bands which have an energy between Elim.
    3) BANDS.get_Ef_from_occ(w_trig=0) returns the value of the last band which an occupation of w_trig.
    4) BANDS.get_bands(bands) returns the projections read in the different PROCAR files (as a CONTRIB object) of the different bands in the "bands" array.
    5) BANDS.get_st(bands) returns a list of different ST objects containing information for plotting the spin texture of the different bands in the "bands" array.
    """    
    
    def __init__(self, path = None, filename = 'PROCAR', filename_outcar = 'OUTCAR', prec_band = 'PROCAR', dlines = 0):
        """
        The class permits to assemble the information of different PROCAR files for different paths. The inputs are :
        - path : the path to the PROCAR and OUTCAR files. If path is an array of different path then BANDS will comute the information of the different paths.
        - filename : name of the PROCAR file. If filename is a string then all the PROCAR files in the path array have to be named as filename. If filename is an array of the same dimension than path then filename_i will coincide with path_i.
        - filename_outcar : name of the OUTCAR file. If filename is a string then all the OUTCAR files in the path array have to be named as filename_outcar. If filename_outcar is an array of the same dimension than path then filename_i will coincide with path_i.
        """
        self.dlines = dlines
        if type(path) == str :
            self.path = [clear_path(path)]
            self.filename = [filename]
            self.filename_outcar = [filename_outcar]
            info = self.get_info_procar()
            self.nkpts = info['nkpts']
            self.nions = info['nions']
            self.nbands = info['nbands']
            self.orbs = info['orbs']
            self.lsorbit = get_lsorbit(self.path[0]+self.filename_outcar[0])
            try :
                self.ispin = get_ISPIN(path)
            except :
                self.ispin = 1
            #self.energy = get_Energy_from_OUTCAR(self.path).reshape((self.nkpts, self.nbands, 3))
            self.energy = self.get_band_energy(prec_band)
            self.kpoints = KPOINTS(self.path[0],  fileformat = 'OUTCAR')
            self.kpoints_dir = KPOINTS(self.path[0],  fileformat = 'OUTCAR', reciprocal = False)
            self.empty = False
        elif type(path) == type(None) :
            self.path = []
            self.filename = []
            self.nkpts = 0
            self.nions = 0
            self.nbands = 0
            self.orbs = []
            self.lsorbit = None
            self.ispin = None
            self.energy = np.array([])
            self.kpoints = KPOINTS(reciprocal = True)
            self.kpoints_dir = KPOINTS(reciprocal = False)
            self.empty = True
        elif type(path) == list :
            new = BANDS(path = None)
            for i, p in enumerate(path) :
                if type(filename) == str :
                    fn = filename
                else :
                    fn = filename[i]
                if type(filename_outcar) == str :
                    fno = filename_outcar
                else :
                    fno = filename_outcar[i]
                tmp = BANDS(p, fn,fno, prec_band)
                new = new + tmp
                
                #for d in tmp1.__dict__.keys() :
                #    new.__dict__[d] = tmp1.__dict__[d]
            for d in new.__dict__.keys() :
                self.__dict__[d] = new.__dict__[d]
    def comp(self, other) :
        C = True
        if self.nions != other.nions :
            C = False
            raise Exception('number of ions are different')
        if self.nbands != other.nbands :
            C = False
            raise Exception('number of bands are different')
        if np.all(self.orbs != other.orbs) :
            raise Exception('orbitals projections are different')
            C = False
        if self.lsorbit != other.lsorbit :
            raise Exception('LSORBIT are different')
            C = False
        if self.ispin != other.ispin :
            raise Exception('ISPIN are different')
            C = False
        return C
        
    
    def __add__(self, other):
        if self.empty == True :
            new = copy.deepcopy(other)
            return new
        if other.empty == True :
            new = copy.deepcopy(self)
            return new
        if self.comp(other) :
            new = copy.deepcopy(self)
            new.path += other.path
            new.filename += other.filename
            new.nkpts += other.nkpts
            new.energy = np.append(self.energy, other.energy, axis = 1)
            new.kpoints += other.kpoints
            new.kpoints_dir += other.kpoints_dir
            return new
    
    def get_band_energy(self, prec_band) :
        """
        Return an nkpts x nbands x 3D np.array containing the band index, the energy and the occupation value of band obtained from the OUTCAR (prec_band = 'OUTCAR') or from the PROCAR file (prec_band = 'PROCAR')
        """
        if prec_band == 'OUTCAR' :
            if self.ispin == 1 :
                E = get_Energy_from_OUTCAR(self.path[0], self.filename_outcar[0]).reshape((1,self.nkpts, self.nbands, 3))
            else :
                E = get_Energy_from_OUTCAR(self.path[0], self.filename_outcar[0]).reshape((2,self.nkpts, self.nbands, 3))
        elif prec_band == 'PROCAR' :
            command = 'grep "^band" %s' %(self.path[0]+self.filename[0])
            f = os.popen(command)
            sent = f.read()
            f.close()
            if self.ispin == 1 :
                E = np.float64(np.array(sent.split()).reshape(1,self.nkpts, self.nbands, 8)[:,:,:,1::3])
            else :
                E = np.float64(np.array(sent.split()).reshape(2,self.nkpts, self.nbands, 8)[:,:,:,1::3])
        return E
    
    def get_bands_limit(self, Elim, Ef = 0,**kwargs):
        """
        Return a 1D np.array giving the bands index whith an energie-Ef included in Elim.
        """
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
            
    
    
    def get_Ef_from_occ(self, w_trig = 0.5):
        """
        Return a float corresponding to the Fermi energy obtained as the last occupied band with an occupation value inferior to w_trig.
        """
        ixs = np.where(self.energy[:,:,:,2] > w_trig)
        return np.max(self.energy[ixs[0],ixs[1],ixs[2],1])
    
    def get_info_procar(self):
        """
        Return a dictionnary containing some informations from the PROCAR files such as :
            - The number of kpoints:  'nkpts' (int)
            - The number of bands:   'nbands' (int)
            - The number of ions:     'nions' (int)
            - The projected orbitals:  'orbs' (1D np.array of string)
        """
        command = "head -2 %s | tail -1" %(self.path[0]+self.filename[0])
        f = os.popen(command)
        sent = f.read()
        f.close()
        data = np.int64(np.array(sent.split())[3::4])
        
        command = "head -8 %s | tail -1" %(self.path[0]+self.filename[0])
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
            a = str(ATOMS[0])
            if a == '0' :
                a = 'tot'
            keys_atoms = '%s' %(a)
        else :
            a = str(ATOMS[0])
            if a == '0' :
                a = 'tot'
            keys_atoms = '(%s' %(a)
            for at in ATOMS[1:] :
                a = str(at)
                if a == '0' :
                    a = 'tot'
                keys_atoms += '|%s' %(a)
            keys_atoms += ')'
            
        if self.ispin == 1 and self.lsorbit == True :
            nlines = 2 + (self.nions+1+self.dlines)*4
            SPINS = ['x','y','z','tot']
        else :
            nlines = 2 + (self.nions+1)
            if self.ispin == 1 :
                SPINS = ['tot']
            else :
                SPINS = ['up', 'dn']
        data_string = ""
        for i,p in enumerate(self.path) :
            command = 'egrep "^band *%s " %s -A %d | egrep "^ *%s "' %(keys_bands, self.path[i]+self.filename[i], nlines, keys_atoms)
            f = os.popen(command)
            data_string += f.read().replace('tot','0')
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
    
    def get_st(self, bands, ATOMS = 0, orbs = 'tot', direct = True) :
        """
        Return an 1D-array of ST objects corresponding to the spin texture of various bands projected on the atoms containing in ATOMS and on orbitals containing on orbs.
        """
        import re
        if type(orbs) == str :
            orbs = re.split('\+|,',orbs)
        if type(ATOMS) == str :
            ATOMS = re.split('\+|,',ATOMS)
        CONTRIBS = self.get_bands(bands)
        STs = []
        for b in range(len(bands)) :
            ix = np.array(self.kpoints.data.drop_duplicates(['kx','ky','kz']).index)
            if direct :
                kpts = np.array(self.kpoints_dir.data[['kx','ky','kz']])[ix]*2*np.pi
            else :
                kpts = np.array(self.kpoints.data[['kx','ky','kz']])[ix]
            sx = CONTRIBS['x'][b][ATOMS][orbs][ix]
            sy = CONTRIBS['y'][b][ATOMS][orbs][ix]
            sz = CONTRIBS['z'][b][ATOMS][orbs][ix]
            S = np.column_stack((sx,sy,sz))
            E = self.energy[0,ix,bands[b]-1,1]
            st = ST(E,S,kpts)
            STs.append(st)
        return STs
    

    
class ST():
    """
    This class gives some function to plot some data of the spin texture such as 
    :
    1) ST.plot   : plot the spin orientation in form of 2D arrows as a function of 2D wavevectors.
    2) ST.plot_S : plot the value of the spin along one direction in color as a function of 2D wavectors.
    3) ST.plot_E : plot the energy in color as a function of 2D wavectors.
    4) ST.ecut   : plot the spin orientation in form of 2D arrows for a given value of the energy as a function of 2D wavevectors.
    """
    def __init__(self, E = None, S = None, K = None,
                 **kwargs) :
        """
        - ST.E    : Energy data (1D array)
        - ST.S    : Spin orientation data (3D array)
        - ST.K    : Wavevectors (3D array)
        - ST.norm : Norm of the spins (1D array)
        """
        
        self.E = E
        self.S = S
        self.K = K
        if type(self.S) == type(None):
            self.norm = None
        else :
            self.norm = np.linalg.norm(self.S, axis = 1)
            
    def plot(self,ix_K = [0,1], ix_S = [0,1,2], plot_sz = True ,N_interp = 0, 
               ax = None, space = 1, 
                   **kwargs):
        """
        1) ST.plot plots the spin orientation in form of 2D arrows as a function of 2D wavevectors. The inputs are :
            - ix_K     : 2D array which gives the 2D wavectors direction as the x and y axis of the graphic (0 = x ; 1 = y ; 2 = z).
            - ix_S     : 2D array which gives the 2D arrows (0 = x ; 1 = y ; 2 = z).
            - plot_sz  : plot the ix_S[2] direction with color if True.
            - N_interp : Size of the interpolation.
            - ax       : pyplot axe to plot.
            - space    : space between 2 plotted points.
            - kwargs   : extra arguments for the ax.quiver function (cmap, norm, color, scale, headwidth, headlength, headaxislength, width, ...). Exemple of kwargs if plot_sz = True : 
            ST.plot(cmap = 'seismic', norm = plt.Normalize(-1,1), scale = 50, headwidth = 3, headlength = 3, headaxislength = 2, width = 0.006)
            Exemple of kwargs if plot_sz = False : 
            ST.plot(color = 'red', scale = 50, headwidth = 3, headlength = 3, headaxislength = 2, width = 0.006)
        """
        if ax == None :
            ax = plt.gca()
        
        
        kx = self.K[:,ix_K[0]]
        ky = self.K[:,ix_K[1]]
        sx = self.S[:,ix_S[0]] / self.norm
        sy = self.S[:,ix_S[1]] / self.norm
        sz = self.S[:,ix_S[2]] / self.norm
        
        if N_interp == 0 :
            if plot_sz : 
                l = ax.quiver(kx[::space],ky[::space], 
                          sx[::space], sy[::space], sz[::space],
                          scale_units='xy',angles = 'xy',**kwargs)
            else :
                l = ax.quiver(kx[::space],ky[::space], 
                          sx[::space], sy[::space],
                          scale_units='xy',angles = 'xy',**kwargs)
            
        else : 
            X = np.linspace(np.min(kx), np.max(ky), N_interp)
            Y = np.linspace(np.min(kx), np.max(ky), N_interp)
            x,y  = np.meshgrid(X, Y)
            
            sx = griddata((kx,ky), sx, (x,y))
            sy = griddata((kx,ky), sy, (x,y))
            sz = griddata((kx,ky), sz, (x,y))
            
            if plot_sz : 
                l = ax.quiver(x[::space,::space],y[::space,::space], 
                          sx[::space,::space], sy[::space,::space],
                          scale_units='xy',angles = 'xy',**kwargs)
            else :
                l = ax.quiver(x[::space,::space],y[::space,::space], 
                          sx[::space,::space], sy[::space,::space], sz[::space,::space],
                          scale_units='xy',angles = 'xy',**kwargs)
        return l 
    
    def plot_S(self,ix_S = 0,ix_K = [0,1], N_interp = 101, 
               ax = None, 
                   **kwargs):
        """
        2) ST.plot_S plots the value of the spin along one direction in color as a function of 2D wavectors. The inputs are :
            - ix_S     : int which gives the plotted spin component (0 = x ; 1 = y ; 2 = z).
            - ix_K     : 2D array which gives the 2D wavectors direction as the x and y axis of the graphic (0 = x ; 1 = y ; 2 = z).
            - N_interp : Size of the interpolation.
            - ax       : pyplot axe to plot.
            - kwargs   : extra arguments for the ax.imshow function (cmap, norm, ...). Exemple of kwargs : 
            ST.plot(cmap = 'seismic', norm = plt.Normalize(-1,1))
        """
        if ax == None :
            ax = plt.gca()
        
        
        kx = self.K[:,ix_K[0]]
        ky = self.K[:,ix_K[1]]
        s = self.S[:,ix_S] / self.norm
        
            
        X = np.linspace(np.min(kx), np.max(ky), N_interp)
        Y = np.linspace(np.min(kx), np.max(ky), N_interp)
        x,y  = np.meshgrid(X, Y)

        s = griddata((kx,ky), s, (x,y))
        l = ax.imshow(s, interpolation = 'nearest', extent = [np.min(x), np.max(x), np.min(y), np.max(y)],**kwargs)
        return l 
            
    def plot_E(self,Ef = 0,ix_K = [0,1], ix_S = [0,1,2],N_interp = 101, 
               ax = None, 
                   **kwargs):
        """
        3) ST.plot_E plots the energy in color as a function of 2D wavectors. The inputs are :
            - Ef   : Fermi energy (default = 0).
            - ix_K     : 2D array which gives the 2D wavectors direction as the x and y axis of the graphic (0 = x ; 1 = y ; 2 = z).
            - ix_S     : 2D array which gives the 2D arrows (0 = x ; 1 = y ; 2 = z).
            - N_interp : Size of the interpolation.
            - ax       : pyplot axe to plot.
            - kwargs   : extra arguments for the ax.imshow function (cmap, norm, ...). Exemple of kwargs : 
            ST.plot_E(Ef = 2.37, cmap = 'gist_rainbow', norm = plt.Normalize(-1,1))
        """
        if ax == None :
            ax = plt.gca()
        
        
        kx = self.K[:,ix_K[0]]
        ky = self.K[:,ix_K[1]]
        sx = self.S[:,ix_S[0]] / self.norm
        sy = self.S[:,ix_S[1]] / self.norm
        sz = self.S[:,ix_S[2]] / self.norm
        
        E = self.E-Ef
            
        X = np.linspace(np.min(kx), np.max(ky), N_interp)
        Y = np.linspace(np.min(kx), np.max(ky), N_interp)
        x,y  = np.meshgrid(X, Y)

        sx = griddata((kx,ky), sx, (x,y))
        sy = griddata((kx,ky), sy, (x,y))
        sz = griddata((kx,ky), sz, (x,y))
        e = griddata((kx,ky), E, (x,y))
        l = ax.imshow(e, interpolation = 'nearest', extent = [np.min(x), np.max(x), np.min(y), np.max(y)],**kwargs)
        return l 
    
    def ecut(self, Ecut, Ef= 0, N_interp = 51, 
             ix_K = [0,1],
             ix_S = [0,1],
             ax = None,
             args_cs = {},
             **kwargs
            ) :
        """
        4) ST.ecut plots the spin orientation in form of 2D arrows for a given value of the energy as a function of 2D wavevectors. The inputs are :
            - Ecut : Cut energy where the spins are plotted.
            - Ef   : Fermi energy (default = 0).
            - N_interp : Size of the interpolation.
            - ix_K     : 2D array which gives the 2D wavectors direction as the x and y axis of the graphic (0 = x ; 1 = y ; 2 = z).
            - ix_S     : 2D array which gives the 2D arrows (0 = x ; 1 = y ; 2 = z).
            - ax       : pyplot axe to plot.
            - args_cs  : Dictionnary for extra arguments of the ax.contour function (colors, linestyle, ...)
            - kwargs   : extra arguments for the ax.quiver function (color, scale, headwidth, headlength, headaxislength, width, ...). Exemple of kwargs :
            ST.ecut(4.2, Ef = 2.37, color = 'red', scale = 50, headwidth = 3, headlength = 3, headaxislength = 2, width = 0.006)
        """
        if 'color' in kwargs.keys() :
            args_cs['colors'] = kwargs.get('color', None)
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
        
        kx = self.K[:,ix_K[0]]
        ky = self.K[:,ix_K[1]]
        sx = self.S[:,ix_S[0]] / self.norm
        sy = self.S[:,ix_S[1]] / self.norm
        
        if ax == None :
            ax = plt.gca()
        X = np.linspace(np.min(kx), np.max(kx), N_interp)
        Y = np.linspace(np.min(ky), np.max(ky), N_interp)
        x,y  = np.meshgrid(X, Y)
        E_grid = griddata((kx,ky), self.E, (x,y))

        CS = ax.contour(x,y,E_grid-Ef,Ecut, **args_cs)
        paths = CS.collections[0].get_paths()
        verts = [xx.vertices for xx in paths]
        points = np.concatenate(verts)
        kx_grid = points[:,0]
        ky_grid = points[:,1]
        sx_grid = griddata((kx,ky), sx, (kx_grid,ky_grid))
        sy_grid = griddata((kx,ky), sy, (kx_grid,ky_grid))

        ax.quiver(kx_grid,ky_grid,sx_grid,sy_grid, **kwargs)
        return kx_grid,ky_grid,sx_grid,sy_grid
    
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
   