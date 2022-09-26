#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 10:42:10 2019

@author: gosteau
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

class HAMILTONIAN():
    class EIG():
        def __init__(self, eigval, eigvec, kx,ky,kz):
            self.eigval = eigval
            self.eigvec = eigvec
            self.n = len(self.eigval)
            self.kx = kx
            self.ky = ky
            self.kz = kz
        
        def __repr__(self):
            string = 'k: %.8f   %.8f   %.8f\n' %(self.kx, self.ky, self.kz)
            for i in range(self.n):
                string += 'E= %4.8f ; v%d: [ ' %(self.eigval[i], i)
                for j in range(self.n-1):
                    string += '{num.real: 10.5f} {num.imag:+10.5f}j  ;  '.format(num=self.eigvec[i][j])
                string += '{num.real: 10.5f} {num.imag:+10.5f}j  ]\n'.format(num=self.eigvec[i][self.n-1])
            return string
        
        def av_val(self,O):
            """
            The function will return an array of the different average value of the Operator for each eigenvectors.
                                    av_val = [av_val_i for each i]
                                    av_val_i =  <psi,i| O |psi,i>
            where |psi,i> is the eigenvector number i ; O is the operator and av_val_i the average value.
            Inputs :
                - O : Operator (matrix)
            """
            tab = []
            for i in range(self.n):
                av_val_i = np.dot(np.conjugate(self.eigvec[i]),np.dot(O, self.eigvec[i]))
                tab.append(av_val_i)
            return np.array(tab)
            
    
    class EIGS():
        def __init__(self,eig):
            self.eig = eig
            self.n = self.eig[0].n
            self.N = len(eig)
            """
            self.eigval = []
            self.eigvec = []
            self.KX = []
            self.KY = []
            self.KZ = []
            for eig in tab_eig:
                self.eigval.append(eig.eigval)
                self.eigvec.append(eig.eigvec)
                self.KX.append(eig.kx)
                self.KY.append(eig.ky)
                self.KZ.append(eig.kz)
            self.eigval = np.array(self.eigval).transpose()
            self.eigvec = np.array(self.eigvec)
            self.KX = np.array(self.KX)
            self.KY = np.array(self.KY)
            self.KZ = np.array(self.KZ)
            """
        def eigval(self):
            return np.array([self.eig[k].eigval for k in range(self.N)]).transpose()
        
        def eigvec(self):
            tmp = np.array([self.eig[k].eigvec for k in range(self.N)])
            return tmp.transpose()
            tmp = np.swapaxes(tmp,0,1)
            return tmp
            return np.swapaxes(tmp,1,2)
        
        def KX(self):
            return np.array([self.eig[k].kx for k in range(self.N)])
        def KY(self):
            return np.array([self.eig[k].ky for k in range(self.N)])
        def KZ(self):
            return np.array([self.eig[k].kz for k in range(self.N)])
        
        def av_val(self, O):
            return np.array([self.eig[k].av_val(O) for k in range(self.N)]).transpose()
    
    def __init__(self, Ham):
        self.Ham = Ham
        
    def resolve(self,kx,ky,kz):
        RES = np.linalg.eig(self.Ham(kx,ky,kz))
        eigval = RES[0]
        eigvec = np.transpose(RES[1])
        #eigvec = RES[1]
        A = np.dot(self.Ham(kx,ky,kz),RES[1])
        B = np.linalg.inv(RES[1])
        diagHam = np.dot(B, A)
        eig = self.EIG(np.sort(eigval), eigvec[np.argsort(eigval)], kx,ky,kz)
        return eig
    
    def resolve_line(self, K: np.array):
        Tab = []
        for k in range(len(K[0])):
            kx = K[0][k]
            ky = K[1][k]
            kz = K[2][k]
            Tab.append(self.resolve(kx,ky,kz))
        return self.EIGS(Tab)
        
    def proj(self,v1,v2, kx,ky,kz):
        return np.dot(np.transpose(v1), np.dot(self.Ham(kx,ky,kz), v2))
    
        
        
    def diagproj(self,v1,v2):
        return np.dot(np.transpose(v1), np.dot(self.diagHam, v2))


