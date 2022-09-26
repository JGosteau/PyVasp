from lib_poscar import *



class POLA_prev():
    def verif(self):
        return True
    
    class LAT():
        def __init__(self):
            self.lattice = None
            self.lattice_name = ''
            self.pion = None
            self.pelc = None
            self.P = None
            self.quanta = None
        
        def __repr__(self):
            conv = self.conv
            unit = self.unit
            string = 'Along %s axis:\n' %(self.lattice_name)
            string += "% 50s" %("Ionic dipole moment: p[ion]=(  ")
            for i in self.pion :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            string += "% 50s" %("Total electronic dipole moment: p[ion]=(  ")
            for i in self.pelc :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            string += "% 50s" %("Total dipole moment: P=(  ")
            for i in self.P :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            string += "% 50s" %("Polarization Quantum: n=(  ")
            for i in self.quanta :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            return string
    
    def get_quanta(self):
        self.cartesian.quanta = np.zeros(3)
        for i in range(3) :
            for j in range(3):
                self.cartesian.quanta[i] += self.direct.lattice[j][i]
    
    def read_pola(self):
        outcar = open(self.path+'OUTCAR','r')
        for line in outcar :
            if "p[ion]" in line :
                self.reciprocal.pion = np.float64(line.split('(')[1].split(')')[0].split())
            elif "p[elc]" in line :
                self.reciprocal.pelc = np.float64(line.split('(')[1].split(')')[0].split())
        outcar.close()
        self.reciprocal.P = self.reciprocal.pion + self.reciprocal.pelc
        
    
    def conversion(self, unit_name = 'electrons Angst'):
        C = 1.6e-19
        A = 1e-10
        V = np.dot(self.direct.lattice[:,0], np.cross(self.direct.lattice[:,1],self.direct.lattice[:,2]))
        superscript_2 = "\u00B2"
        superscript_minus = "\u207B"
        if unit_name == 'SI' or unit_name == 'C.m-2':
            conv = C*A/(V*A**3)
            unit_name = 'C.m'+superscript_minus+superscript_2
        elif unit_name == 'uC.cm-2' :
            conv = C*1e6*(A*1e2)/(V*(A*1e2)**3)
            unit_name = '\u03BC'+'C.cm'+superscript_minus+superscript_2
        else :
            conv = 1
            unit_name = 'electrons Angst'
        
        self.conv = conv
        self.unit = unit_name
        
        self.reciprocal.conv = conv
        self.reciprocal.unit = unit_name
        
        self.direct.conv = conv
        self.direct.unit = unit_name
        
        self.cartesian.conv = conv
        self.cartesian.unit = unit_name
            
    
    def __init__(self, path):
        self.path = path
        self.verif()
        info = INFO_OUTCAR(self.path)
        
        self.reciprocal = self.LAT()
        #self.reciprocal.lattice = np.transpose(info.rec_lattice/info.length_rec_vectors)
        self.reciprocal.lattice = np.transpose(info.rec_lattice)/info.length_rec_vectors
        self.reciprocal.lattice_name = 'reciprocal lattice'
        self.direct = self.LAT()
        self.direct.lattice_name = 'direct lattice'
        self.direct.lattice = np.transpose(info.direct_lattice)
        self.cartesian = self.LAT()
        self.cartesian.lattice_name = 'cartesian'
        
        self.read_pola()
        self.get_quanta()
        
        self.cartesian.pion = np.dot(self.reciprocal.lattice, self.reciprocal.pion)
        self.cartesian.pelc = np.dot(self.reciprocal.lattice, self.reciprocal.pelc)
        self.cartesian.P = np.dot(self.reciprocal.lattice, self.reciprocal.P)
        
        #d_lattice_normed = np.transpose(info.direct_lattice/info.length_direct_vectors)
        d_lattice_normed = np.transpose(info.direct_lattice)/info.length_direct_vectors
        self.direct.pion = np.dot(np.linalg.inv(d_lattice_normed), self.cartesian.pion)
        self.direct.pelc = np.dot(np.linalg.inv(d_lattice_normed), self.cartesian.pelc)
        self.direct.P = np.dot(np.linalg.inv(d_lattice_normed), self.cartesian.P)
        self.direct.quanta = np.dot(np.linalg.inv(d_lattice_normed), self.cartesian.quanta)
        
        self.reciprocal.quanta = np.dot(np.linalg.inv(self.reciprocal.lattice), self.cartesian.quanta)
        self.conversion()
    
    def __repr__(self):
        string = "%s\n" %(str(self.cartesian))
        string += "%s\n" %(str(self.direct))
        string += "%s\n" %(str(self.reciprocal))
        return string
    
class POLA():
    def verif(self):
        return True
    
    class LAT():
        def __init__(self):
            self.lattice = None
            self.lattice_name = ''
            self.pion = None
            self.pelc = None
            self.P = None
            self.quanta = None
        
        def __repr__(self):
            conv = self.conv
            unit = self.unit
            string = 'Along %s axis:\n' %(self.lattice_name)
            string += "% 50s" %("Ionic dipole moment: p[ion]=(  ")
            for i in self.pion :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            string += "% 50s" %("Total electronic dipole moment: p[ion]=(  ")
            for i in self.pelc :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            string += "% 50s" %("Total dipole moment: P=(  ")
            for i in self.P :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            string += "% 50s" %("Polarization Quantum: n=(  ")
            for i in self.quanta :
                string += " % 5.5f " %(i*conv) 
            string += " ) %s\n" %(unit)
                
            return string
    
    def get_quanta(self):
        self.cartesian.quanta = np.zeros(3)
        for i in range(3) :
            for j in range(3):
                self.cartesian.quanta[i] += self.direct.lattice[j][i]
    
    def read_pola(self):
        outcar = open(self.path+'OUTCAR','r')
        for line in outcar :
            if "p[ion]" in line :
                self.cartesian.pion = np.float64(line.split('(')[1].split(')')[0].split())
            elif "p[elc]" in line :
                self.cartesian.pelc = np.float64(line.split('(')[1].split(')')[0].split())
        outcar.close()
        print(self.cartesian.pion)
        print(self.cartesian.pelc)
        self.cartesian.P = self.cartesian.pion + self.cartesian.pelc
        
    
    def conversion(self, unit_name = 'electrons Angst'):
        C = 1.6e-19
        A = 1e-10
        V = np.dot(self.direct.lattice[:,0], np.cross(self.direct.lattice[:,1],self.direct.lattice[:,2]))
        superscript_2 = "\u00B2"
        superscript_minus = "\u207B"
        if unit_name == 'SI' or unit_name == 'C.m-2':
            conv = C*A/(V*A**3)
            unit_name = 'C.m'+superscript_minus+superscript_2
        elif unit_name == 'uC.cm-2' :
            conv = C*1e6*(A*1e2)/(V*(A*1e2)**3)
            unit_name = '\u03BC'+'C.cm'+superscript_minus+superscript_2
        else :
            conv = 1
            unit_name = 'electrons Angst'
        
        self.conv = conv
        self.unit = unit_name
        
        self.reciprocal.conv = conv
        self.reciprocal.unit = unit_name
        
        self.direct.conv = conv
        self.direct.unit = unit_name
        
        self.cartesian.conv = conv
        self.cartesian.unit = unit_name
            
    
    def to_cart(self, V, B):
        #P = np.zeros(3)
        #for i in range(3) :
        #    v = V[i]
        #    u = B[i]
        #    P += v*u
        #return P
        P = np.dot(np.transpose(B),V)
        return P
    
    def to_base(self, V,B) :
        P = np.dot(np.linalg.inv(np.transpose(B)), V)
        return P
            
    
    def __init__(self, path):
        self.path = path
        self.verif()
        info = INFO_OUTCAR(self.path)
        
        self.reciprocal = self.LAT()
        #self.reciprocal.lattice = np.transpose(info.rec_lattice/info.length_rec_vectors)
        self.reciprocal.lattice = info.rec_lattice#/info.length_rec_vectors
        #for i in range(3) :
        #    self.reciprocal.lattice[i] /= info.length_rec_vectors[i]
        self.reciprocal.lattice_name = 'reciprocal lattice'
        self.direct = self.LAT()
        self.direct.lattice_name = 'direct lattice'
        self.direct.lattice = info.direct_lattice
        self.cartesian = self.LAT()
        self.cartesian.lattice_name = 'cartesian'
        
        self.read_pola()
        self.direct.quanta = info.length_direct_vectors
        self.direct.normed = copy.deepcopy(self.direct.lattice)
        for i in range(3) :
            self.direct.normed[i] /= info.length_direct_vectors[i]
        self.reciprocal.normed = copy.deepcopy(self.reciprocal.lattice)
        for i in range(3) :
            self.reciprocal.normed[i] /= info.length_rec_vectors[i]
        
        
        #self.get_quanta()
        #self.cartesian.pion = self.to_cart(self.reciprocal.pion, self.reciprocal.lattice) 
        #self.cartesian.pelc = self.to_cart(self.reciprocal.pelc, self.reciprocal.lattice) 
        #self.cartesian.P = self.to_cart(self.reciprocal.P, self.reciprocal.lattice) 
        self.cartesian.quanta = self.to_cart(self.direct.quanta,  self.direct.normed)
                
        self.direct.pion = self.to_base(self.cartesian.pion, self.direct.normed)
        self.direct.pelc = self.to_base(self.cartesian.pelc, self.direct.normed)
        self.direct.P = self.to_base(self.cartesian.P, self.direct.normed)
        #self.direct.quanta = self.to_base(self.cartesian.quanta,  self.direct.lattice/info.length_direct_vectors)
        
        self.reciprocal.pion = self.to_base(self.cartesian.pion, self.reciprocal.normed)
        self.reciprocal.pelc = self.to_base(self.cartesian.pelc, self.reciprocal.normed)
        self.reciprocal.P = self.to_base(self.cartesian.P, self.reciprocal.normed)
        self.reciprocal.quanta = self.to_base(self.cartesian.quanta,  self.reciprocal.lattice)
        self.conversion()
    
    def __repr__(self):
        string = "%s\n" %(str(self.cartesian))
        string += "%s\n" %(str(self.direct))
        string += "%s\n" %(str(self.reciprocal))
        return string

class BEC():
    """
    Class which stocks the different values of the Born Effective Charges $\bm{Z^*_i}$ calculated by VASP using LEPSILON = .TRUE. 
    The inner parameters are :
        - path:  path to directory where filename is present
        - filename: name of the OUTCAR file (default = "OUTCAR")
    """
    
    
    def get_data(self):
        n = grep('BORN', self.path + self.filename)[-1][-1]
        tmp = POSCAR(self.path)
        #self.n_atoms = len(grep(' ion   ', self.path + self.filename)[-1])
        self.n_atoms = len(tmp)
        npos = 14
        string = get_row(n+2, self.path + self.filename, self.n_atoms*4)
        data = {}
        for at in range(self.n_atoms):
            data[at+1] = np.float64(np.array(string.split())[2+npos*at:2+npos*at+12].reshape(3,4)[:,1:])
        return data
    
    def __init__(self, path, filename = 'OUTCAR', DATAS = None):
        if type(DATAS) == type(None) :
            self.filename = filename
            self.path = path
            self.data = self.get_data()
        else :
            self.filename = 'Datas'
            self.path = None
            data = np.array(list(DATAS.values()))
            self.data = {}
            for at in range(len(data)):
                self.data[at+1] = data[at]
        
    def __getitem__(self, ix):
        return self.data[ix]
    
    def __repr__(self):
        string = ""
        for at in self.data.keys() : 
            string += "ion %3d:\n" %(at)
            string += str(self[at]) + '\n'
        return string
    
    def conv(self, unit):
        if unit == 'C.A':
            unit = 1.6e-19
        elif unit == 'C.m-2':
            unit = 1.6e-19/((1e-10)**2)
        elif unit == 'uC.cm-2':
            unit = 1.6e-19*1e6/((1e-8)**2)
        else :
            unit = 1
        return unit
    
    def delta_P_at(self, atom, poscar_i, poscar_f, unit = 'uC.cm-2') :
        """
        Permit to calculate the variation of polarization from 1 atom between position 1 to position 2.
        \delta \bm{P}(n) = \frac{e}{V} \sum_i \bm{Z^*}(n) \cdot \delta\bm{r}(n)        
        where n is the atom, e the electric charge and V the volume of the cell
        """
        unit = self.conv(unit)
        f1 = unit/poscar_i.V
        f2 = unit/poscar_f.V
        x2,y2,z2 = np.array(poscar_f.data.loc[atom, ['x','y','z']])
        pos_f = x2*poscar_f.lattice[0] + y2*poscar_f.lattice[1] + z2*poscar_f.lattice[2]
        
        x1,y1,z1 = np.array(poscar_i.data.loc[atom, ['x','y','z']])
        pos_i = x1*poscar_i.lattice[0] + y1*poscar_i.lattice[1] + z1*poscar_i.lattice[2]
        
        x = ((x2-x1)+0.5)%1-0.5
        y = ((y2-y1)+0.5)%1-0.5
        z = ((z2-z1)+0.5)%1-0.5
        
        pos = x*poscar_f.lattice[0] + y*poscar_f.lattice[1] + z*poscar_f.lattice[2] 
        #dr = np.float64(pos_f*f2- pos_i*f1)        
        dr = np.float64(pos*f2)
        return np.dot(self.data[atom], np.transpose(dr))
        
    def delta_P(self, poscar_i, poscar_f, unit = 'uC.cm-2') :
        P = np.zeros(3)
        for at in self.data.keys() :
            P += self.delta_P_at(at, poscar_i, poscar_f)
        return P
    
    
    