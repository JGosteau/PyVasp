from lib_gen import *


class wannier_Exchange_int():
    def read(self, ini):
        with open(self.file, 'r') as f :
            string = ""
            n = 0
            while n < ini :
                f.readline()
                n +=1
            for line in f.readlines() :
                string += line
        data = np.float64(string.split())
        #data.reshape((self.n_points, self.n_wannier**2, self.n_format_int))
        data = data.reshape((self.n_points * self.n_wannier**2,self.n_format_int ))
        return data
    
    def get_positions(self) :
        pos = grep('unit_cell_cart',self.path+self.seedname+'.win')[1]
        self.lattice = np.float64(get_row(pos[0]+1, self.path+self.seedname+'.win', a = 2).split()).reshape((3,3))
        a,b,c = self.lattice
        V = np.dot(a, np.cross(b,c))
        u = np.cross(b,c)/V
        v = np.cross(c,a)/V
        w = np.cross(a,b)/V
        self.rec_lat = np.array([u,v,w])*2*np.pi
    
    def __init__(self, path, seedname = "wannier90", suffix = "_hr.dat", n_format_deg = 15, n_format_int = 7):
        self.path = clear_path(path)
        self.seedname = seedname
        self.file = self.path + seedname + suffix
        self.n_wannier  = int(get_row(2,self.file).split()[0])
        self.n_points  = int(get_row(3,self.file).split()[0])
        self.n_format_int = n_format_int
        n_lines = int(np.ceil(self.n_points/n_format_deg))
        self.deg = np.int64(get_row(4,self.file, n_lines-1).split())
        if len(self.deg) != self.n_points :
            self.deg = np.int64(get_row(4,self.file, n_lines).split())
        self.data = self.read(4+n_lines-1)
        self.posx = np.unique(self.data[:,0])
        self.posy = np.unique(self.data[:,1])
        self.posz = np.unique(self.data[:,2])
        self.kpts = []
        
        self.IX = []
        for n1 in range(1,self.n_wannier+1):
            tmp = []
            for n2 in range(1,self.n_wannier+1):
                tmp.append(self.get_coord_orb_b(n1,n2))
            self.IX.append(tmp)
        self.IX = np.array(self.IX)
        self.get_positions()
        
        
    def get_coord_atom(self,x,y,z) :
        ix = np.where((self.data[:,0] == x) * (self.data[:,1] == y) * (self.data[:,2] == z))
        return self.data[ix]
    
    def get_coord_orb_b(self,n1,n2) :
        ix = np.where((self.data[:,3] == n1) * (self.data[:,4] == n2))[0]
        return ix
    
    def get_coord_orb(self,n1,n2) :
        ix = np.where((self.data[:,3] == n1) * (self.data[:,4] == n2))
        return self.data[ix]
    
    def get_coord_orb(self,n1,n2) :
        ix = self.IX[n1-1,n2-1]
        return self.data[ix]
    
    def __getitem__(self, ix) :
        return self.get_coord_atom(*ix)
    
    def Sum_all(self) :
        REAL = np.zeros((self.n_wannier, self.n_wannier))
        IM = np.zeros((self.n_wannier, self.n_wannier))
        for n1 in range(1,self.n_wannier + 1):
            for n2 in range(1,self.n_wannier + 1):
                MAT = self.get_coord_orb(n1,n2)[:,[5,6]]
                REAL[n1-1, n2-1] = np.sum(MAT[:,0]/self.deg)
                IM[n1-1, n2-1] = np.sum(MAT[:,1]/self.deg)
        return REAL, IM
    
    def read_kpts_wannier(self, filename1 = None, filename2 = None, path = None) :
        if path == None :
            path = self.path
        if filename1 == None :
            f1 = self.seedname + "_band.kpt"
            filename1 = path+f1
        if filename2 == None :
            f2 = self.seedname + "_band.dat"
            filename2 = path+f2
        string = ""
        with open(filename1, 'r') as f:
            nkpt = int(f.readline().split()[0])
            for line in f.readlines() :
                string += line
        norm = np.loadtxt(filename2, max_rows = nkpt)       
        return np.column_stack((np.float64(string.split()).reshape((nkpt, 4)), norm[:,0]))
    
    
    def get_Wmn_k(self, Kpts):
        n = self.n_wannier
        nkpts = len(Kpts)
        EIG = [] #np.zeros((nkpts,n))
        R_kpts = []
        I_kpts = []
        for ik,K in enumerate(Kpts) :
            REAL = np.zeros((n, n))
            IM = np.zeros((n, n))     
            for n1 in range(n):
                for n2 in range(n):
                    re,im = self.FT(K[[0,1,2]], n1+1,n2+1)
                    REAL[n1,n2] = re
                    IM[n1,n2] = im
            R_kpts.append(REAL)
            I_kpts.append(IM)
        return np.array(R_kpts), np.array(I_kpts)
    
    
    
    def get_EIGS_from_KPTS(self, Kpts, R_kpts, I_kpts, used_col = None):
        if type(used_col) == type(None) :
            used_col = list(range(self.n_wannier))
        n = len(used_col)
        nkpts = len(Kpts)
        EIG = [] #np.zeros((nkpts,n))
        for ik,K in enumerate(Kpts) :
            REAL = R_kpts[ik][used_col][:,used_col]
            IM = I_kpts[ik][used_col][:,used_col]
            eig = np.linalg.eigvals(REAL+1j*IM)
            EIG.append(eig)
        return np.array(EIG)
    
    def get_EIGVECTS_from_KPTS(self, Kpts, R_kpts, I_kpts, used_col = None):
        if type(used_col) == type(None) :
            used_col = list(range(self.n_wannier))
        n = len(used_col)
        nkpts = len(Kpts)
        EIG = [] #np.zeros((nkpts,n))
        VECT = []
        for ik,K in enumerate(Kpts) :
            REAL = R_kpts[ik][used_col][:,used_col]
            IM = I_kpts[ik][used_col][:,used_col]
            eig, vect = np.linalg.eig(REAL+1j*IM)
            ix = np.argsort(eig)
            #print(eig.shape)
            vect = np.transpose(vect)
            EIG.append(eig[ix])
            VECT.append(vect[ix])
        EIG = np.array(EIG)
        VECT = np.array(VECT)
        #ix = np.argsort(EIG)
        #EIG = EIG[ix]
        #VECT = VECT[ix]
        return EIG, VECT
            
            
    def FT(self,K, n1,n2) :
        tmp = self.get_coord_orb(n1,n2)
        R = tmp[:,[0,1,2]]
        f = (tmp[:,-2]+1j*tmp[:,-1])/self.deg
        S = np.sum(f*np.exp(-1j*2*np.pi*np.dot(R,K)))
        return np.real(S),np.imag(S)
    
    def read_bands(self, filename = "wannier90_band.dat"):
        with open(self.path+filename, 'r') as f :
            string = f.read()
        tmp = string.split()
        data = np.float64(tmp).reshape(self.n_wannier, int(len(tmp)/self.n_wannier/2),2)
        return data
    
    def Fourier_transform(self,K) :
        REAL = np.zeros((self.n_wannier, self.n_wannier))
        IM = np.zeros((self.n_wannier, self.n_wannier))        
        for n1 in range(1,self.n_wannier + 1):
            for n2 in range(1,self.n_wannier + 1):
                f = (self.get_coord_orb(n1,n2)[:,-2]+1j*self.get_coord_orb(n1,n2)[:,-1])/self.deg
                R = self.get_coord_orb(n1,n2)[:,[0,1,2]]
                S = np.sum(f*np.exp(-1j*np.dot(R,K)))
                REAL[n1-1,n2-1] = np.real(S)
                IM[n1-1,n2-1] = np.imag(S)
        return REAL, IM
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    