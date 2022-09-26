import sympy as sp
import numpy as np
sp.init_printing(use_unicode=True, use_latex='mathjax')

class HAMILTONIAN(): 
    def get_coeff(self, coeff = None):
        self.coeff = None
        free_symbols = list(self.H.free_symbols)
        for i in range(3):
            try :
                free_symbols.remove(self.K[i])
            except :
                None
        from copy import deepcopy
        ref_symbols = deepcopy(free_symbols)
        if type(coeff) == type(None):
            self.coeff = free_symbols
        else :
            for c in coeff :
                try :
                    free_symbols.remove(c)
                except :
                    raise Exception('Error: %s not in %s' %(str(c), str(ref_symbols)))
            if len(free_symbols) != 0 :
                raise Exception('Error: %s not in %s' %(str(free_symbols), str(coeff)))
            else :
                self.coeff = coeff
            
    
    def __init__(self, Matrix, coeff = None):
        self.H = Matrix
        self.K = sp.symbols("k_x k_y k_z", real = True)
        self.get_coeff(coeff = coeff)
        self.func = sp.lambdify((self.coeff,self.K), self.H)
        self.det_HE()
        self.eigen_mul = None
        self.eigenvals = None
        self.eigenvects = None
    
    def get_eigen(self):
        RES = self.H.eigenvects()
        self.eigen_mul = np.array(RES)[:,1]
        self.eigenvals= np.array(RES)[:,0]
        self.eigenvects = np.array(RES)[:,2]
    
    
    def det_HE(self):
        E = sp.symbols('E', real = True)
        self.det = sp.det(self.H - E*sp.eye(self.H.shape[0])).expand()
        tmp_func = sp.lambdify((self.coeff,self.K,E), self.det)
        self.det_func = lambda coeff, K, E : np.array([tmp_func(coeff,K[i], E[i]) for i in range(len(E))])
    
    def leastsq(self, K,E, ini_val = None, fixed_val = None):
        from scipy.optimize import leastsq
        from copy import deepcopy
        def get_coeff(new_coeff, fixed_val):
            coeff = np.array([])
            i=0
            count = 0
            for j in fixed_val.keys() :
                coeff = np.append(coeff,new_coeff[i:j-count])
                val = fixed_val[j]
                coeff = np.append(coeff, val)
                i = j-count
                count += 1
            coeff = np.append(coeff,new_coeff[i:])
            return coeff
        
        if type(fixed_val) != type(None) :
            Fvals = {}
            new_ini_val = np.delete(ini_val,fixed_val)
            for i in np.sort(fixed_val) :
                Fvals[i] = ini_val[i]
            
            det = lambda new_coeff, K, E : self.det_func(get_coeff(new_coeff, Fvals),K,E)
            res = leastsq(det, new_ini_val, args = (K,E), full_output = 1)
            return res
            
        else :
            det = self.det_func
            res = leastsq(det, ini_val, args = (K,E), full_output = 1)
            return res
        """
        if len(E) == self.H.shape[0] :
            #det = lambda coeff, K, E : np.real(np.sum([self.det_func(coeff,K,E[i]) for i in range(len(E))], axis = 1))
            det = lambda coeff, K, E : np.real(np.append(self.det_func(coeff,K,E[0]), [ self.det_func(coeff,K,E[i]) for i in range(1,len(E))]))
        else :
            if len(E) == 1 :
                E = E[0]
        """
        
    
    def __repr__(self):
        return self.H