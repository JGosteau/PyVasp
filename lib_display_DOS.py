#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 16:01:15 2019

@author: gosteau
"""


from lib_DOS_v2 import *



def display_LDOS(path, atoms = None, Ef = None, doscar = None, orbitals = None, ax = None, spin = None, Elim = None, scale = 1.,  **kwargs):
    if doscar == None:
        doscar = DOSCAR(path, atoms, orbitals = orbitals, Elim = Elim)
    Ef = doscar.Ef
    if orbitals == None:
        orbitals = doscar.used_orbitals
    if ax == None :
        ax = plt.subplot()
    if spin == None :
        if doscar.spin == 2 :
            spin = '_up'
        elif doscar.lsorbit :
            spin = '_sz'
        else :
            spin = ''
    RES = np.zeros(len(doscar.energy))
    for l in orbitals :
        RES += doscar.data[l+spin]
    ax.plot(doscar.energy-Ef, RES*scale, **kwargs)
    return doscar 
    
def display_LDOS_filled(path, atoms = None, Ef = None, doscar = None, orbitals = None, ax = None, spin = None, Elim = None,scale = 1.,  **kwargs):
    if doscar == None:
        doscar = DOSCAR(path, atoms, orbitals = orbitals, Elim = Elim)
    Ef = doscar.Ef
    if orbitals == None:
        orbitals = doscar.used_orbitals
    if ax == None :
        ax = plt.subplot()
    if spin == None :
        if doscar.spin == 2 :
            spin = '_up'
        elif doscar.lsorbit :
            spin = '_sz'
        else :
            spin = ''
    RES = np.zeros(len(doscar.energy))
    for l in orbitals :
        RES += doscar.data[l+spin]
    ax.fill_between(doscar.energy['Energy']-Ef, RES*scale, **kwargs)
    return doscar



   
def multi_disp(path ,poscar = None,plane = None, Ef = None,name = 'dos',orbitals = ['tot'], plane_range = None, n_plane = 5, Elim = None, ylim = None, figsize = (6,12), spins = {'none' : 1}, count_scale = 1, print_label = True, label_plane = False, IF = 0, labelpad = 20, SbH = False, seuil_SbH = 0.1):
    if poscar == None :
        poscar = POSCAR(path, Cond_pos_at=True)
    if plane == None :
        plane = poscar.get_plane()
    if plane_range == None:
        plane_range = range(len(plane))
    FIG = []
    count = 0
    columns = ['name']
    for sp in spins.keys():
            if sp == 'none':
                if dos.spin == 1:
                    sp = ''
                else :
                    sp = '_up'
            for l in orbitals :
                Ec = "Ec(%s)" %(l+sp)
                Ev = "Ev(%s)" %(l+sp)
                Eg = "Eg(%s)" %(l+sp)
                columns += [Ev, Ec, Eg]
    SBH = pd.DataFrame(index = plane_range, columns = columns)
    for p in plane_range :
        #print "plane ",p," : " ,
        if count%n_plane == 0:
            fig = plt.figure(name+" "+str(p)+" --> "+str(p+n_plane), figsize= figsize)
            FIG.append(fig)
            
        ax = plt.subplot(n_plane,1,count%n_plane+1)
        
        dos = DOSCAR(path, atoms = plane[p], orbitals = orbitals, Elim = Elim, Ef = Ef)
        for sp in spins.keys():
            sc = spins[sp]
            if sp == 'none':
                if dos.spin == 1:
                    sp = ''
                else :
                    sp = '_up'
            display_LDOS_filled(path, plane[p],orbitals=['tot'],doscar = dos,ax = ax,Elim = Elim,color = 'grey', alpha = 0.1, spin = sp, scale = sc)
            display_LDOS(path, plane[p],orbitals=['tot'],doscar = dos,ax = ax,Elim = Elim,color = 'black', linestyle = '--', spin = sp, scale = sc)
        if ylim != None: 
            plt.ylim(*ylim)
        if Elim != None:
            plt.xlim(*Elim)
        plt.vlines(0, *plt.ylim(), linestyles=':', color = 'black')
        A = poscar.data.loc[plane[p]]
        unique, counts = np.unique(np.array(A['At']), return_counts=True)
        for at in unique :
            color = COLORS[at]
            atoms = np.array(A['n_at'].loc[A['At'] == at])
            #print at+" :", atoms, ";",
        for sp in spins.keys():
            sc = spins[sp]
            if sp == 'none':
                if dos.spin == 1:
                    sp = ''
                else :
                    sp = '_up'
            display_LDOS(path, plane[p],orbitals=['tot'],doscar = dos,ax = ax,Elim = Elim,color = color, linestyle = '-', spin = sp, scale = sc)
        #print " "
        string = ""
        if 'dxy' in orbitals :
            display_LDOS(path, plane[p],orbitals=['dxy'],doscar = dos,ax = ax,Elim = Elim,color = 'green', linestyle = '-')
            display_LDOS(path, plane[p],orbitals=['dxz','dyz'],doscar = dos,ax = ax,Elim = Elim,color = 'purple', linestyle = '-')
        
        def order_atom(atoms, count, order) :
            tmp = []
            for at in atoms :
                i = 0
                for comp in order :
                    if at == comp :
                        tmp.append(i)
                        break
                    i += 1
            return atoms[np.array(tmp).argsort()], count[np.array(tmp).argsort()]
        def int_comp(a):
            if int(a) == a:
                return int(a)
            else :
                return a
        
        ordered, count_ordered = order_atom(unique, counts, poscar.order)
        for i in range(len(ordered)):
            string += ordered[i]
            val = int_comp(count_ordered[i]/float(count_scale))
            if val != 1 :
                string += str(val)
        def get_label_plane(p):
            if label_plane == True :
                new = (IF-p+1)/2
                if new > 0:
                    return 'IF+'+str(new)
                elif new == 0 :
                    return 'IF-'+str(new)
                else :
                    return 'IF'+str(new)
            else :
                return str(p)
            
        print( string + " " + get_label_plane(p))
        SBH.loc[p]['name'] = string + " " + get_label_plane(p)
        if SbH :
            
            for sp in spins.keys():
                if sp == 'none':
                    if dos.spin == 1:
                        sp = ''
                    else :
                        sp = '_up'
                
                for l in orbitals :
                    orb = l+sp
                    Ec, Ev = dos.SbH(orb, seuil_SbH)
                    Ec_name = "Ec(%s)" %(orb)
                    Ev_name = "Ev(%s)" %(orb)
                    Eg_name = "Eg(%s)" %(orb)
                    SBH.loc[p][Ev_name] = Ev
                    SBH.loc[p][Ec_name] = Ec
                    SBH.loc[p][Eg_name] = Ec-Ev
                    print ("Ev(%s) = %.2f |" %(orb,Ev),end = '')
                    print ("Ec(%s) = %.2f ; " %(orb,Ec))
            print( "--------------")
        if p%2 == 0:
            color = 'white'
        else :
            color = 'grey'
            
            
        box = dict(facecolor=color, pad=5, alpha=0.2)
        if print_label :
            plt.ylabel(string + " " + get_label_plane(p), va  = 'center', ha = 'center', rotation = 0, labelpad = labelpad, bbox = box)
        plt.tick_params(axis="x",which = 'major',direction="in", length=3.0,width = 1.5)
        plt.tick_params(axis="x",which = 'minor',direction="in", length=3.0,width = 1.0)
        plt.tick_params(axis="x",which = 'major',top = True,direction="in", length=3.0,width = 1.5)
        plt.tick_params(axis="x",which = 'minor',top = True, direction="in", length=3.0,width = 1.0)
        plt.tick_params(axis='y',left = False)
        plt.yticks(visible= False)
        plt.xticks(visible= False)
        count += 1
        
        if count %n_plane == 0:
            plt.xticks(visible= True)
            fig.tight_layout()
            plt.subplots_adjust(wspace=0, hspace=0)
        
    return FIG, SBH





#path = "/home/gosteau/Documents/These/1ere année/LAOSTO/NP_IF/new_relax/"
