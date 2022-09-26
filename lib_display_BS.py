#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 09:53:10 2019

@author: gosteau
"""
from lib_class_bs import *

def get_contrib_max(A) :
    B = []
    K = []
    for a in A.keys() :
        K.append(a)
        B.append(A[a])
    ix = np.argmax(B)
    return K[ix]

def get_ener_contrib(band,kpt,Ef=0): 
    D = {}
    E = band.energy[kpt]
    for i in range(len(band.orbs[:-1])) :
        orb = band.orbs[i]
        val = band.tokpt(kpt)[0][i]
        D[orb] = val
    print(E-Ef,get_contrib_max(D))
    return E-Ef,get_contrib_max(D)





def Band_structure(path, Elim, Ef, kpt, knorm = None, picker = 1):
    bands = BANDS(path)
    band_account = bands.get_bands_in_Elim(Elim, Ef)


def Band_Structure_interp(path, Elim, Ef, kpt, knorm = None,
                          orb = 'tot', atom = 'tot', spin = '',
                          interp_kind = None,n_interp = 1000, 
                          BANDS = None, Disp_BANDS = None, 
                          cmap = plt.get_cmap('seismic'), norm = None, picker = 1 ):
    """
    Permet de tracer la stucture de bande et faire des projections sur les atomes et orbitales.
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi 
    - kpt : index des k-points à afficher
    - knorm : norme des directions affichées (absisse). Par défault prend l'index des points kpt
    
    - orb = 'tot' : orbitale projetée ('tot', 'px', 'px+py', ...)
    - at = 'tot' : atome projeté ('tot', 1, '1', '1+2', ...)
    - spin = '' : projection sur les spins '' / '_sx' / '_sy' / '_sz'
    
    - interp_kind : type d'interpolation si None --> pas d'interpolation ('cubic' / 'linear' / None)
    - n_interp : nombre de point pour l'interpolation.

    - BANDS : dictionnaire répertoriant les objets CONTRIB des différentes bandes déjà calculés.
    None par défault (utile lorsqu'on trace plusieurs fois les mêmes bandes).
    - Disp_BANDS : Liste des bandes à afficher. Définit par Elim par défault.
    
    - cmap : color_map utilisée.
    - norm : norme pour la projection.
    - picker : Tolérance pour selectionner un point.
    
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Elim = [Emin, Emax]
        
        Ef = get_Ef(path, fileformat = 'PROCAR')
        kpoints = KPOINTS(path)
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = 1, fay = 1, faz = 1) # (knorm en unité 2pi/a)
        
        poscar = POSCAR(path)
        fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz) # (knorm en unité A-1)
        
        Band_Structure_interp(path, Elim, Ef, kpt, knorm)
    """
    
    def string_spin(spin):
        if 'x' in spin :
            return '_sx'
        elif 'y' in spin :
            return '_sy'
        elif 'z' in spin :
            return '_sz'
        elif 'up' in spin :
            return '_up'
        elif 'dn' in spin :
            return '_dn'
        else :
            return ''
    
    if norm == None and ( 'x' in spin or 'y' in spin or 'z' in spin) :
        norm = mpl.colors.Normalize(-1,1)    
    elif norm == None :
        norm = mpl.colors.Normalize(0,1)  
    
    if BANDS == None :
        BANDS = {}
    if type(Disp_BANDS) == type(None) :
        Disp_BANDS = get_bands_limit(path, Elim, Ef)        
    for b in Disp_BANDS:
        if b not in BANDS.keys() :
            bands = CONTRIB(path, b, string_spin(spin))
            if spin =='_up' or spin == '_dn' :
                BANDS[str(b)+spin] = bands
            else :
                BANDS[str(b)] = bands
        else : 
            bands = BANDS[b]
        k = range(len(kpt))
        new_knorm = knorm
        if interp_kind != None :
            new_bands = inter.interp1d(k,bands.energy[kpt], kind = interp_kind)
            new_orb = inter.interp1d(k,bands[atom][orb][kpt], kind = interp_kind)
            new_norm = inter.interp1d(k,knorm, kind = interp_kind)
            k_new = np.linspace(k[0],k[-1],n_interp)
            band_out = new_bands(k_new)
            orb_out = new_orb(k_new)
            k = k_new
            new_knorm = new_norm(k)
            
        else :
            band_out = bands.energy[kpt]
            orb_out = bands[atom][orb][kpt]
        
        
        if type(knorm) == type(None) :
            index = k
            index_def = range(len(kpt))
        else :
            index = new_knorm
            index_def = knorm
            
        points = np.array([index,band_out-Ef]).T.reshape(-1,1,2)
        segments = np.concatenate([points[:-1], points[1:]], axis = 1)
        lc = LineCollection(segments,cmap=cmap,norm = norm)
        lc.set_array(orb_out)
        plt.gca().add_collection(lc)
        plt.plot(index_def,bands.energy[kpt]-Ef,c = 'black', alpha = 0, linewidth = 1, label = str(b)+bands.spin, picker = picker, zorder = 10)
        plt.ylim(*Elim)
    return BANDS


def Band_Structure_interp_point_orbs(path, Elim, Ef, kpt, knorm = None,
                          ORBS = {}, COLORS_ORBS = {}, spin = '',
                          interp_kind = None,n_interp = 1000, 
                          Disp_BANDS = None,plot_label = True,
                          ax = None,
                          picker = 1 , step = 1, size = 100, alpha = 0.5,
                          **kwargs):
    """
    Permet de tracer la stucture de bande et faire des projections sur les atomes et orbitales.
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi 
    - kpt : index des k-points à afficher
    - knorm : norme des directions affichées (absisse). Par défault prend l'index des points kpt
    
    - orb = 'tot' : orbitale projetée ('tot', 'px', 'px+py', ...)
    - at = 'tot' : atome projeté ('tot', 1, '1', '1+2', ...)
    - spin = '' : projection sur les spins '' / '_sx' / '_sy' / '_sz'
    
    - interp_kind : type d'interpolation si None --> pas d'interpolation ('cubic' / 'linear' / None)
    - n_interp : nombre de point pour l'interpolation.

    - BANDS : dictionnaire répertoriant les objets CONTRIB des différentes bandes déjà calculés.
    None par défault (utile lorsqu'on trace plusieurs fois les mêmes bandes).
    - Disp_BANDS : Liste des bandes à afficher. Définit par Elim par défault.
    
    - cmap : color_map utilisée.
    - norm : norme pour la projection.
    - picker : Tolérance pour selectionner un point.
    
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Elim = [Emin, Emax]
        
        Ef = get_Ef(path, fileformat = 'PROCAR')
        kpoints = KPOINTS(path)
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = 1, fay = 1, faz = 1) # (knorm en unité 2pi/a)
        
        poscar = POSCAR(path)
        fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz) # (knorm en unité A-1)
        
        Band_Structure_interp(path, Elim, Ef, kpt, knorm)
    """
    if ax == None :
        ax = plt.gca()
    def string_spin(spin):
        if 'x' in spin :
            return '_sx'
        elif 'y' in spin :
            return '_sy'
        elif 'z' in spin :
            return '_sz'
        elif 'up' in spin :
            return '_up'
        elif 'dn' in spin :
            return '_dn'
        else :
            return ''
    
    if type(Disp_BANDS) == type(None) :
        Disp_BANDS = get_bands_limit(path, Elim, Ef)        
    BANDS = {}
    handles = []
    for b in Disp_BANDS:
        bands = CONTRIB(path, b, string_spin(spin))
        BANDS[b] = bands
        k = range(len(kpt))
        new_knorm = knorm
        band_out = bands.energy[kpt]
        
        
        if type(knorm) == type(None) :
            index = k
            index_def = range(len(kpt))
        else :
            index = new_knorm
            index_def = knorm
            
        #plt.plot(index_def, bands.energy[kpt]-Ef, c= 'black', label = str(b)+bands.spin, picker = picker)
        #ax.plot(index_def, bands.energy[kpt]-Ef, c= color, picker = picker, alpha = alpha_band, **kwargs)
        
        def get_atoms_and_orbs(l):
            def get_at(tmp_at):
                tmp = tmp_at.split(',')
                atoms = []
                for ats in tmp :
                    if '-' in ats :
                        tmp2 = ats.split('-')
                        ini = int(tmp2[0])
                        end = int(tmp2[-1])+1
                        step = 1
                        atoms += list(range(ini,end, step))
                    else :
                        atoms += [int(ats)]
                return atoms
            tmp = l.split('_')
            if len(tmp) == 2 :
                tmp_at = tmp[0]#.split(',')
                atoms = get_at(tmp_at)
                orb = tmp[1]
            else :
                if tmp[0][0] in '1234567890' :
                    atoms = get_at(tmp[0])
                    orb = 'tot'
                else :
                    atoms = 0
                    orb = tmp[0]
            return atoms, orb
        
        for l in COLORS_ORBS.keys():
            atoms, orb = get_atoms_and_orbs(l)
            s = bands[atoms][orb][kpt]*size
            marker = 'o'
            c = COLORS_ORBS[l]
            if l not in handles and plot_label == True:
                label_plot = l
                ax.plot([np.nan],[np.nan], color = c, markersize = 10, label = label_plot, alpha = alpha, marker = marker, linewidth = 0)
                handles.append(l)
            else :
                label_plot = None
            ax.scatter(index_def[::step], bands.energy[kpt][::step]-Ef, c = c, s = s[::step], marker = marker, alpha = alpha, **kwargs)
        
        ax.set_ylim(*Elim)
    return BANDS

def Band_Structure(path, Elim, Ef, kpt, knorm= None,
                  interp_kind = None,n_interp = 1000, spin = '',
                  Disp_BANDS = None, ax = None,color = 'black', **kwargs) :
    if ax == None :
        ax = plt.gca()
    def string_spin(spin):
        if 'x' in spin :
            return '_sx'
        elif 'y' in spin :
            return '_sy'
        elif 'z' in spin :
            return '_sz'
        elif 'up' in spin :
            return '_up'
        elif 'dn' in spin :
            return '_dn'
        else :
            return ''
    if type(Disp_BANDS) == type(None) :
        Disp_BANDS = get_bands_limit(path, Elim, Ef)  
    bands,occ  = get_band_info(path, spin = string_spin(spin))
    for b in Disp_BANDS:
        band = bands[b-1]
        k = range(len(kpt))
        new_knorm = knorm
        
        if type(knorm) == type(None) :
            index = k
            index_def = range(len(kpt))
        else :
            index = new_knorm
            index_def = knorm
        ax.plot(index_def, band[kpt]-Ef, c= color, **kwargs)
        ax.set_ylim(*Elim)

def Band_Structure_interp_point_spins(path, Elim, Ef, kpt, knorm = None,
                          spin = 'sx', state = 'up', spin_threshold = 0,
                          interp_kind = None,n_interp = 1000, 
                          Disp_BANDS = None,plot_label = True,
                          ax = None,
                          picker = 1 , step = 1, size = 100, alpha = 0.5, color = 'red', marker = 'o',
                          **kwargs):
    if ax == None :
        ax = plt.gca()
    def string_spin(spin):
        if 'x' in spin :
            return '_sx'
        elif 'y' in spin :
            return '_sy'
        elif 'z' in spin :
            return '_sz'
        elif 'up' in spin :
            return '_up'
        elif 'dn' in spin :
            return '_dn'
        else :
            return ''
    if type(Disp_BANDS) == type(None) :
        Disp_BANDS = get_bands_limit(path, Elim, Ef)  
    handles = []
    for b in Disp_BANDS:
        bands = CONTRIB(path, b, string_spin(spin))
        k = range(len(kpt))
        new_knorm = knorm

        if type(knorm) == type(None) :
            index = k
            index_def = range(len(kpt))
        else :
            index = new_knorm
            index_def = knorm
        tmp = bands['tot']['tot'][kpt]
        s = np.full(len(tmp),np.nan)
        if state == 'up' :
            ix = np.where(tmp > np.abs(spin_threshold))[0]
            s[ix] = tmp[ix]*size
        else :
            ix = np.where(tmp < -np.abs(spin_threshold))[0]
            s[ix] = -tmp[ix]*size
        #print(b, bands.energy[kpt], tmp>0, tmp<0)
        ax.scatter(index_def[::step], bands.energy[kpt][::step]-Ef, c = color, s = s[::step], marker = marker, alpha = alpha)

        ax.set_ylim(*Elim)





def Band_Structure_interp_point_orbs_dE(path, Elim, Ef, B1,B2,kpt, knorm = None,
                          ORBS = {}, COLORS_ORBS = {}, spin = '',
                          interp_kind = None,n_interp = 1000, 
                          Disp_BANDS = None, alpha_band = 1,plot_label = True,
                          ax = None,
                          picker = 1 , step = 1, size = 100, alpha = 0.5, color = 'black',
                          **kwargs):
    """
    Permet de tracer la stucture de bande et faire des projections sur les atomes et orbitales.
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi 
    - kpt : index des k-points à afficher
    - knorm : norme des directions affichées (absisse). Par défault prend l'index des points kpt
    
    - orb = 'tot' : orbitale projetée ('tot', 'px', 'px+py', ...)
    - at = 'tot' : atome projeté ('tot', 1, '1', '1+2', ...)
    - spin = '' : projection sur les spins '' / '_sx' / '_sy' / '_sz'
    
    - interp_kind : type d'interpolation si None --> pas d'interpolation ('cubic' / 'linear' / None)
    - n_interp : nombre de point pour l'interpolation.

    - BANDS : dictionnaire répertoriant les objets CONTRIB des différentes bandes déjà calculés.
    None par défault (utile lorsqu'on trace plusieurs fois les mêmes bandes).
    - Disp_BANDS : Liste des bandes à afficher. Définit par Elim par défault.
    
    - cmap : color_map utilisée.
    - norm : norme pour la projection.
    - picker : Tolérance pour selectionner un point.
    
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Elim = [Emin, Emax]
        
        Ef = get_Ef(path, fileformat = 'PROCAR')
        kpoints = KPOINTS(path)
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = 1, fay = 1, faz = 1) # (knorm en unité 2pi/a)
        
        poscar = POSCAR(path)
        fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz) # (knorm en unité A-1)
        
        Band_Structure_interp(path, Elim, Ef, kpt, knorm)
    """
    if ax == None :
        ax = plt.gca()
    def string_spin(spin):
        if 'x' in spin :
            return '_sx'
        elif 'y' in spin :
            return '_sy'
        elif 'z' in spin :
            return '_sz'
        elif 'up' in spin :
            return '_up'
        elif 'dn' in spin :
            return '_dn'
        else :
            return ''
    
    if type(Disp_BANDS) == type(None) :
        Disp_BANDS = get_bands_limit(path, Elim, Ef)        
    BANDS = {}
    handles = []
    for i in range(len(B1)):
        b1 = B1[i]
        b2 = B2[i]
        bands_1 = CONTRIB(path, b1, string_spin(spin))
        bands_2 = CONTRIB(path, b2, string_spin(spin))
        BANDS[b1] = bands_1
        BANDS[b2] = bands_2
        k = range(len(kpt))
        new_knorm = knorm
        band_out_1 = bands_1.energy[kpt]
        band_out_2 = bands_2.energy[kpt]
        band_out = band_out_2 - band_out_1
        
        if type(knorm) == type(None) :
            index = k
            index_def = range(len(kpt))
        else :
            index = new_knorm
            index_def = knorm
            
        #plt.plot(index_def, bands.energy[kpt]-Ef, c= 'black', label = str(b)+bands.spin, picker = picker)
        ax.plot(index_def, band_out, c= color, picker = picker, alpha = alpha_band, **kwargs)
        def get_atoms_and_orbs(l):
            def get_at(tmp_at):
                tmp = tmp_at.split(',')
                atoms = []
                for ats in tmp :
                    if '-' in ats :
                        tmp2 = ats.split('-')
                        ini = int(tmp2[0])
                        end = int(tmp2[-1])+1
                        step = 1
                        atoms += list(range(ini,end, step))
                    else :
                        atoms += [int(ats)]
                return atoms
            tmp = l.split('_')
            if len(tmp) == 2 :
                tmp_at = tmp[0]#.split(',')
                atoms = get_at(tmp_at)
                orb = tmp[1]
            else :
                if tmp[0][0] in '1234567890' :
                    atoms = get_at(tmp[0])
                    orb = 'tot'
                else :
                    atoms = 0
                    orb = tmp[0]
            return atoms, orb
        
        for l in COLORS_ORBS.keys():
            atoms, orb = get_atoms_and_orbs(l)
            print(atoms, orb)
            s = bands_1[atoms][orb][kpt]*size
            marker = 'o'
            c = COLORS_ORBS[l]
            if l not in handles and plot_label == True:
                label_plot = l
                ax.plot([np.nan],[np.nan], color = c, markersize = 10, label = label_plot, alpha = alpha, marker = marker, linewidth = 0)
                handles.append(l)
            else :
                label_plot = None
            ax.scatter(index_def[::step], band_out[::step], c = c, s = s[::step], marker = marker, alpha = alpha)
        
        #ax.set_ylim(*Elim)
    return BANDS

#TAB_pick = []

def plot_BS_contrib(path, Elim, Ef = None, HS_lines = 'auto',
                    kpoints = None, poscar = None, kpt = None, knorm = None,
                    ispin=None, cmap_ispin = None,
                    fig = None, ax = None, figname = 'BS', onpick_bool = True,
                    verbose = False,
                    **args):
    """
    Permet de tracer la stucture de bande et faire des projections sur les atomes et orbitales.
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi. Prend le dernier niveau occupé par défault.
    - HS_lines : Direction de Haute-symétrie, doit être de la forme : 
    HS_lines = [[hs1,hs2], [hs3,hs4], ...] où [hs1,hs2] représente la direction de l'objet HS : hs1 à hs2 (voir help(HS))
    Par défault : 'auto' cherche automatiquement les directions de Hautes Symétries.
    
    - kpoints : objet KPOINTS repertoriant la grille de k-points.
    Par défault prend les kpoints données par kpts_tmp dans le dossier path/BANDS
    généré par la fonction create_bands(path)
    - poscar : objet POSCAR repertoriant la position des atomes.
    Par défault prend la position des atomes donnés dans le OUTCAR du dossier path.
    
    - ispin : permet de tracer les 2 spins sur un même graphe si ISPIN = 2.
    Par défault, le programme cherche si il y a deux spins et les affichent en même temps.
    si ispin = 1, la fonction ne tracera qu'un seul spin definit par le tag spin.
    - cmap_ispin : couleur des diférents spin. doit être de la forme :
    cmap_ispin = {'up' : color_map_up, 'dn' : color_map_dn}
    Par défault les deux spins ont la même couleur.
    
    - fig : objet pyplot.figure
    - ax : objet pyplot.axes
    - figname : donne le nom de la figure si fig == None
    - onpick_bool : définit si les tracer sont sélectionnable (True) où non (False).
    
    - verbose : affiche des informations supplémentaires.
    - **args : arguments de la fonction Band_Structure_interp(**args).
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Elim = [Emin, Emax]
        fig, ax, BANDS = plot_BS_contrib(path, Elim)
    """   
    
    
    
    if ispin == None :
        if len(np.where(np.array(os.listdir(path+'BANDS/')) == 'band_info_tmp_up')[0]) == 1 :
            ispin = 2
        else :
            ispin = 1
    if Ef == None :
        Ef = get_Ef(path, fileformat = 'PROCAR')
    if type(kpoints) == type(None) :
        kpoints = KPOINTS(path)
    if type(poscar) == type(None) :
        poscar = POSCAR(path)
    fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
    if type(HS_lines) == str :
        if HS_lines == 'auto' :
            HS_lines = list(kpoints.auto_lines())
            if verbose :
                print("Found %d HS directions : " %(len(HS_lines)))
                for L in HS_lines :
                    print("   %s -> %s" %(L[0], L[1]))
    elif verbose :
        print("%d HS directions : " %(len(HS_lines)))
        for L in HS_lines :
            print("   %s -> %s" %(L[0], L[1]))
    if type(kpt) == type(None) or type(knorm) == type(None) :
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz)
    if fig == None : 
        fig = plt.figure(figname)
        ax = None
    if ax == None :
        ax = plt.subplot()
        
    if ispin == 2 :
        if cmap_ispin != None :
            try :
                cmap_up = cmap_ispin['_up']
            except :
                cmap_up = cmap_ispin['up']
            try :
                cmap_dn = cmap_ispin['_dn']
            except :
                cmap_dn = cmap_ispin['dn']
            args['cmap'] = cmap_up
        args['spin'] = '_up'
        BANDS = Band_Structure_interp(path, Elim, Ef, kpt, knorm, **args)
        if cmap_ispin != None :
            args['cmap'] = cmap_dn
        args['spin'] = '_dn'
        args['BANDS'] = BANDS
        BANDS = Band_Structure_interp(path, Elim, Ef, kpt, knorm,**args)
    else : 
        BANDS = Band_Structure_interp(path, Elim, Ef, kpt, knorm, **args)
    
    HSs = {}
    try :
        for i in HS_lines :
            if i[0].name not in HSs :
                HSs[i[0].name] =  i[0]
    except :
        HSs = HS_points
    ticks, label = new_ticks(kpoints,kpt,knorm, HSs = HSs)
    plt.xticks(knorm[np.int0(ticks)],label)
    ax.set_xticks(knorm[np.int0(ticks)])
    ax.set_xticklabels(label)
    ax.vlines(knorm[np.int0(ticks)],*plt.ylim())
    ax.set_xlim(knorm[0], knorm[-1])
    ax.hlines(0,*plt.xlim(), color = 'black', linestyle = ':')

    def onpick(event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ldata = thisline.get_label()
        ind = event.ind
        k_ind = ind
        try : 
            band = np.int0(ldata)
        except :
            band = np.int0(ldata.split('_')[0])
        band = ldata
        for i in range(len(k_ind)) :
            index = kpt[k_ind[i]]
            supp_string = " norm %2.5f " %(xdata[i])
            supp_string = None
            BANDS[band].print_tokpt(index, kpoints = kpoints,poscar = poscar, Ef = Ef, supp_string = supp_string)

    if onpick_bool :
        fig.canvas.mpl_connect('pick_event', onpick) 
    
    
    return fig, ax, BANDS

def plot_BS_contrib_point(path, Elim, Ef = None, HS_lines = 'auto',
                    kpoints = None, poscar = None, kpt = None, knorm = None,
                    ispin=None, cmap_ispin = None,
                    fig = None, ax = None, figname = 'BS', onpick_bool = True,
                    verbose = False,
                    **args):
    """
    Permet de tracer la stucture de bande et faire des projections sur les atomes et orbitales.
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi. Prend le dernier niveau occupé par défault.
    - HS_lines : Direction de Haute-symétrie, doit être de la forme : 
    HS_lines = [[hs1,hs2], [hs3,hs4], ...] où [hs1,hs2] représente la direction de l'objet HS : hs1 à hs2 (voir help(HS))
    Par défault : 'auto' cherche automatiquement les directions de Hautes Symétries.
    
    - kpoints : objet KPOINTS repertoriant la grille de k-points.
    Par défault prend les kpoints données par kpts_tmp dans le dossier path/BANDS
    généré par la fonction create_bands(path)
    - poscar : objet POSCAR repertoriant la position des atomes.
    Par défault prend la position des atomes donnés dans le OUTCAR du dossier path.
    
    - ispin : permet de tracer les 2 spins sur un même graphe si ISPIN = 2.
    Par défault, le programme cherche si il y a deux spins et les affichent en même temps.
    si ispin = 1, la fonction ne tracera qu'un seul spin definit par le tag spin.
    - cmap_ispin : couleur des diférents spin. doit être de la forme :
    cmap_ispin = {'up' : color_map_up, 'dn' : color_map_dn}
    Par défault les deux spins ont la même couleur.
    
    - fig : objet pyplot.figure
    - ax : objet pyplot.axes
    - figname : donne le nom de la figure si fig == None
    - onpick_bool : définit si les tracer sont sélectionnable (True) où non (False).
    
    - verbose : affiche des informations supplémentaires.
    - **args : arguments de la fonction Band_Structure_interp(**args).
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Elim = [Emin, Emax]
        fig, ax, BANDS = plot_BS_contrib(path, Elim)
    """   
    
    
    
    if ispin == None :
        if len(np.where(np.array(os.listdir(path+'BANDS/')) == 'band_info_tmp_up')[0]) == 1 :
            ispin = 2
        else :
            ispin = 1
    if Ef == None :
        Ef = get_Ef(path, fileformat = 'PROCAR')
    if type(kpoints) == type(None) :
        kpoints = KPOINTS(path)
    if type(poscar) == type(None) :
        poscar = POSCAR(path)
    fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
    if type(HS_lines) == str :
        if HS_lines == 'auto' :
            HS_lines = list(kpoints.auto_lines())
            if verbose :
                print("Found %d HS directions : " %(len(HS_lines)))
                for L in HS_lines :
                    print("   %s -> %s" %(L[0], L[1]))
    elif verbose :
        print("%d HS directions : " %(len(HS_lines)))
        for L in HS_lines :
            print("   %s -> %s" %(L[0], L[1]))
    if type(kpt) == type(None) or type(knorm) == type(None) :
        kpt, knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz)
    if fig == None : 
        fig = plt.figure(figname)
        ax = None
    if ax == None :
        ax = plt.subplot()
    print(args.keys(),Ef)
    if ispin == 2 :
        args['spin'] = '_up'
        BANDS = Band_Structure_interp_point_orbs(path, Elim, Ef, kpt, knorm, **args)
        args['spin'] = '_dn'
        BANDS = Band_Structure_interp_point_orbs(path, Elim, Ef, kpt, knorm,**args)
    else : 
        BANDS = Band_Structure_interp_point_orbs(path, Elim, Ef, kpt, knorm, **args)
    print(BANDS.keys())
    HSs = {}
    try :
        for i in HS_lines :
            if i[0].name not in HSs :
                HSs[i[0].name] =  i[0]
    except :
        HSs = HS_points
    ticks, label = new_ticks(kpoints,kpt,knorm, HSs = HSs)
    plt.xticks(knorm[np.int0(ticks)],label)
    ax.set_xticks(knorm[np.int0(ticks)])
    ax.set_xticklabels(label)
    ax.vlines(knorm[np.int0(ticks)],*plt.ylim())
    ax.set_xlim(knorm[0], knorm[-1])
    ax.hlines(0,*plt.xlim(), color = 'black', linestyle = ':')

    def onpick(event):
        thisline = event.artist
        xdata = thisline.get_xdata()
        ydata = thisline.get_ydata()
        ldata = thisline.get_label()
        ind = event.ind
        k_ind = ind
        try : 
            band = np.int0(ldata)
        except :
            band = np.int0(ldata.split('_')[0])
        band = ldata
        for i in range(len(k_ind)) :
            index = kpt[k_ind[i]]
            supp_string = " norm %2.5f " %(xdata[i])
            supp_string = None
            BANDS[band].print_tokpt(index, kpoints = kpoints,poscar = poscar, Ef = Ef, supp_string = supp_string)

    if onpick_bool :
        fig.canvas.mpl_connect('pick_event', onpick) 
    
    
    return fig, ax, BANDS

def routine_BS_d_2(path, Elim, Ef = None, HS_lines = 'auto',
                 ORBS = {'dxy' : 'o', 'dyz' : 'o', 'dxz' : 'o', 'dz2' : 'o', 'x2-y2' : 'o'}, COLORS_ORBS = {'dxy' : 'red', 'dyz' : 'lawngreen', 'dxz' : 'deeppink', 'dz2' : 'blue', 'x2-y2' : 'cyan'}, 
                 step = 2, alpha = 0.15,
                 fig = None, ax = None, figname = None, legend = False,bool_tight_layout = False, kpoints = None, poscar = None,
                 **args) :
    """
    Permet de tracer la stucture de bande et faire des projections sur les orbitales.
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi. Prend le dernier niveau occupé par défault.
    
    - seuil : seuil de tolérance pour afficher la projection sur une orbitale.
    
    - fig : objet pyplot.figure
    - ax : objet pyplot.axes
    - figname : donne le nom de la figure si fig == None
    - legend : Affiche (True) ou non (False) la légende
    
    - **args : arguments de la fonction plot_BS_contrib(**args).
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Elim = [Emin, Emax]
        fig, ax, BANDS = routine_BS_d(path, Elim)
    """       
    

    if figname == None: 
        figname = 'BS orb decomp'
    if fig == None :
        fig = plt.figure(figname)
    if ax == None :
        ax = plt.subplot()
        ax.set_ylabel('E-Ef (eV)')
        bool_tight_layout = True
    if Ef == None :
        Ef = get_Ef(path, fileformat = 'PROCAR')
    if type(kpoints) == type(None):
        kpoints = KPOINTS(path, fileformat = 'PROCAR')
    if type(HS_lines) == str :
        if HS_lines == 'auto' :
            HS_lines = list(kpoints.auto_lines())
    if type(poscar) == type(None):
        poscar = POSCAR(path)
    fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
    kpt, knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz)
    args['Disp_BANDS'] = get_bands_limit(path, Elim, Ef)
    args['fig'] = fig
    args['ax'] = ax
    onpick_bool = args.get('onpick_bool', True)
    try :
        args.pop('onpick_bool')
    except :
        None
    args['ORBS'] = ORBS
    args['COLORS_ORBS'] = COLORS_ORBS
    args['step'] = step
    args['alpha'] = alpha
    BANDS = plot_BS_contrib_point(path,Elim, Ef,kpt = kpt, knorm = knorm,HS_lines = HS_lines,**args)
    #custom_lines = [Line2D([0], [0], color=cmap_n(0), label = 'tot')]
    
    if legend :
        fig.legend(handles=custom_lines, loc = 'upper center', ncol = len(ORBS.keys())+1, mode = 'expand')   
    if bool_tight_layout :
        fig.tight_layout()
    return fig,BANDS


def routine_BS_d(path, Elim, Ef = None, HS_lines = 'auto',
                 ORBS = {'dxy' : cmap_r, 'dyz' : cmap_g, 'dxz' : cmap_p, 'dz2' : cmap_b, 'x2-y2' : cmap_c, 'p' : cmap_o},
                 seuil = 0.5,
                 fig = None, ax = None, figname = None, legend = False,bool_tight_layout = False,
                 **args) :
    """
    Permet de tracer la stucture de bande et faire des projections sur les orbitales.
    - path : "chemin vers le dossier où 'BANDS' est present"
    - Elim = [Emin,Emax] : gamme d'énergie à tracer 
    - Ef : niveau de Fermi. Prend le dernier niveau occupé par défault.
    
    - seuil : seuil de tolérance pour afficher la projection sur une orbitale.
    
    - fig : objet pyplot.figure
    - ax : objet pyplot.axes
    - figname : donne le nom de la figure si fig == None
    - legend : Affiche (True) ou non (False) la légende
    
    - **args : arguments de la fonction plot_BS_contrib(**args).
    
    Utilisation :
        path = "chemin vers le dossier où 'BANDS' est present"
        Elim = [Emin, Emax]
        fig, ax, BANDS = routine_BS_d(path, Elim)
    """       
    

    if figname == None: 
        figname = 'BS orb decomp'
    if fig == None :
        fig = plt.figure(figname)
    if ax == None :
        ax = plt.subplot()
        ax.set_ylabel('E-Ef (eV)')
        bool_tight_layout = True
    norm = Normalize(seuil-0.1,seuil)
    if Ef == None :
        Ef = get_Ef(path, fileformat = 'PROCAR')
    kpoints = KPOINTS(path, fileformat = 'PROCAR')
    if type(HS_lines) == str :
        if HS_lines == 'auto' :
            HS_lines = list(kpoints.auto_lines())
    poscar = POSCAR(path)
    fax, fay, faz = np.array(poscar.diag_lattice)/(2*np.pi)
    kpt, knorm = kpoints.get_dirs(HS_lines, fax = fax, fay = fay, faz = faz)
    args['Disp_BANDS'] = get_bands_limit(path, Elim, Ef)
    args['fig'] = fig
    args['ax'] = ax
    onpick_bool = args.get('onpick_bool', True)
    try :
        args.pop('onpick_bool')
    except :
        None
    BANDS = plot_BS_contrib(path,Elim, Ef=Ef, cmap = cmap_n,onpick_bool = onpick_bool,kpt = kpt, knorm = knorm,HS_lines = HS_lines,**args)[-1]
    args['verbose'] = False
    custom_lines = [Line2D([0], [0], color=cmap_n(0), label = 'tot')]
    for orb in ORBS.keys():
        cmap = ORBS[orb]
        custom_lines.append(Line2D([0], [0], color=cmap(0), label = orb))
        if orb == 'p':
            orb = 'px+py+pz'
        BANDS = plot_BS_contrib(path,Elim, Ef=Ef, orb = orb,cmap = cmap,norm = norm,onpick_bool=False, picker = 0,BANDS = BANDS,kpt = kpt, knorm = knorm,HS_lines = HS_lines,**args)[-1]
    if legend :
        fig.legend(handles=custom_lines, loc = 'upper center', ncol = len(ORBS.keys())+1, mode = 'expand')   
    if bool_tight_layout :
        fig.tight_layout()
    return fig,BANDS

# TEST 
#Hlines = ['G','X','M','G']
path = '/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/constrain/a388_c416/SOC/BS/'
path_PTO = '/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/constrain/a388_c416/SOC/BS/ZRAZ/'
#Hlines = [['G','X'],['X','M'],['Z','R'],['A','Z']]
"""
path = "/home/gosteau/Documents/TMP_remi/PTO/Spin_text/test1/constrain/a388_c416/BS/"
fig1 = plt.figure('up')
fig1, ax, bands = plot_BS_contrib(path, [-3,5], cmap = cmap_n, interp_kind = None, spin = '_up', onpick_bool=True, fig = fig1)


fig2 = plt.figure('dn')
fig2, ax, bands = plot_BS_contrib(path, [-3,5], cmap = cmap_n, interp_kind = None, spin = '_dn', onpick_bool=True, fig = fig2)
#ax = plt.twinx()
#ax = plt.subplot(1,2,2)
#fig,ax, bands = plot_BS_contrib(path, [-3,5], cmap = cmap_r, interp_kind = None, spin = '_d', fig = fig, ax = ax,onpick_bool = False)



path_wannier = path+'wannier/'
Ef= get_Ef(path,fileformat = 'PROCAR')
band, kpt = get_bands_from_wannier(path_wannier, spin = 'up')
for b in band :
    fig1.axes[0].plot(kpt, b-Ef, color = 'red')
    
band, kpt = get_bands_from_wannier(path_wannier, spin = 'dn')
for b in band :
    fig2.axes[0].plot(kpt, b-Ef, color = 'red')
"""
