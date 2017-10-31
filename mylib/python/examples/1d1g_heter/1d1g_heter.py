# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:22:51 2016

@author: segonpin
"""
# we add the root directory for the python library
# if we do not have that directory in our path
import sys
sys.path.insert(0, '../..')
# we use numpy for the data
#import numpy as np
#
from prob.runforest import RunForest
from prob.settings import setmethod, setalgebra
#from prob.vtkparser1d import plot_data

def test():
    """
    >>> PROB_NAME = "1d1g_heter"
    >>> forest = RunForest(prob_name=PROB_NAME)
    >>> forest_refl = RunForest(prob_name=PROB_NAME+"_refl")
    >>> METHOD = setmethod(mtype='transport', dsa=False, quad=['GaussLegendre', 96])
    >>> ALGEBRA = setalgebra(matfree=True)
    >>> FE_SETTINGS= dict([('degree', 2), ('n_ref', 2)])
    >>> forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    >>> forest_refl.set_settings(method=METHOD, fe_settings=FE_SETTINGS)
    >>> LATT = dict()
    >>> LATT[0] = dict([('type', 'grid'), \
                    ('nnodes', [7]), \
                    ('components', [0, 1, 0, 1, 0, 1, 0]), \
                    ('length', [2.7, 2.4, 2.7, 2.4, 2.7, 2.4, 2.7])])
    >>> CORE = dict([('composed', True), \
                 ('nnodes', [1]), \
                 ('length', [[18.0]]), \
                 ('components', [0]), \
                 ('boundary', [0, 0])]) 
    >>> forest.set_geometry(dim=1, core=CORE, lattices=LATT)
    >>> LATT = dict()
    >>> LATT[0] = dict([('type', 'grid'), \
                    ('nnodes', [4]), \
                    ('components', [1, 0, 1, 0]), \
                    ('length', [1.2, 2.7, 2.4, 2.7])])
    >>> CORE = dict([('composed', True), \
                 ('nnodes', [1]), \
                 ('length', [[9.0]]), \
                 ('components', [0]), \
                 ('boundary', [2, 0])])
    >>> forest_refl.set_geometry(dim=1, core=CORE, lattices=LATT)
    >>> MIX = dict()
    >>> MIX[0] = dict([('name', 'moderator'), \
                   ('SigmaT', [0.37037000]), \
                   ('Chi',    [1.00000000]), \
                   ('SigmaS', [0.33400000]), \
                   ('NuSigF', [0.00000000])])
    >>> MIX[1] = dict([('name', 'fuel'), \
                   ('SigmaT', 0.41666700), \
                   ('Chi',    1.00000000), \
                   ('SigmaS', 0.33400000), \
                   ('NuSigF', 0.17800000)])
    >>> forest.set_materials(ngroups=1, mix=MIX)
    >>> forest_refl.set_materials(ngroups=1, mix=MIX)
    >>> forest.run() # doctest: +ELLIPSIS
    Running: ...
    >>> forest_refl.run() # doctest: +ELLIPSIS
    Running: ...
    >>> keff = forest.get_keff()
    >>> keff_refl = forest_refl.get_keff()
    >>> abs(keff-1.16222540417) < 1.e-6
    True
    >>> abs(keff_refl-1.16222594896) < 1.e-6
    True
    >>> abs(keff - keff_refl)< 1.e-6
    True
    """
    # -------------------------------------------------------------------------
    TESTDIR = "1d1g_heter/"
    PROB_NAME = "1d1g_heter"
    DESCRIPTION = (
        "Seven Regions Heterogeneous 1D Slab \n \n"
        "The reference result is keff = 1.16224, and it has \n"
        "been obtained from ONEDANT, using an angular \n"
        "quadrature order of S96, with 500 fine mesh cells \n"
        "in each region and a convergence criterion of \n"
        "tol = 1.e-6, as it was reported in: \n \n"
        "M. Capilla,  CF. Talavera,  D. Ginestar,  G. Verdu,  \n"
        "\"A nodal collocation method for the calculation of \n"
        "the lambda modes of the P_L equations'; \n"
        "Annals of Nuclear Energy, 32,  (2005), 1825-1853.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    forest_refl = RunForest(prob_name=PROB_NAME+"_refl", description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport', dsa=False, quad=['GaussLegendre', 96])
    #METHOD = setmethod(mtype='transport', quad=['GaussLegendre', 96])
    #METHOD = setmethod(mtype='diffusion')
    #ALGEBRA = setalgebra(matfree=True)
    #ALGEBRA = setalgebra(matfree=True, eig=["PI",1.e-7,1000])
    #ALGEBRA = setalgebra(matfree=True, form=,'Standard',eig=["Krylov:1:3",1.e-7,1000], inner=["Krylov",1.e-9,1000,1])
    ALGEBRA = setalgebra(matfree=True, form='Standard',eig=["Krylov:1:3",1.e-7,1000], mg=["krylov",1.e-7,100], inner=["Krylov",1.e-9,1000])
    #ALGEBRA = setalgebra(matfree=True, eig=["PI",1.e-7,1000,'Standard'], inner=["Krylov",1.e-9,1000,1])
    FE_SETTINGS= dict([('degree', 2), ('n_ref', 2)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    forest_refl.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    # Geometry ----------------------------------------------------------------
    LATT = dict()
    lenr = 2.7 #1.0/0.371
    lenf = 2.4 #1.0/0.415
    LATT[0] = dict([('type', 'grid'),
                    ('nnodes', [7]),
                    ('components', [0, 1, 0, 1, 0, 1, 0]),
                    ('length', [lenr, lenf, lenr, lenf, lenr, lenf, lenr])])
    CORE = dict([('composed', True),
                 ('nnodes', [1]),
                 ('length', [[18.0]]),
                 ('components', [0]),
                 ('boundary', [0, 0])])
    forest.set_geometry(dim=1, core=CORE, lattices=LATT)
    
    # Geometry ----------------------------------------------------------------
    LATT = dict()
    lenr = 2.7 #1.0/0.371
    lenf = 2.4 #1.0/0.415
    LATT[0] = dict([('type', 'grid'),
                    ('nnodes', [4]),
                    ('components', [0, 1, 0, 1]),
                    ('length', [lenr, lenf, lenr, lenf/2])])
    CORE = dict([('composed', True),
                 ('nnodes', [1]),
                 ('length', [[9.0]]),
                 ('components', [0]),
                 ('boundary', [0, 2])])
    forest_refl.set_geometry(dim=1, core=CORE, lattices=LATT)
    
    # Materials ---------------------------------------------------------------
    MIX = dict()
    MIX[0] = dict([('name', 'moderator'),
                   ('SigmaT', [0.37037000]),
                   ('Chi',    [1.00000000]),
                   ('SigmaS', [0.33400000]),
                   ('NuSigF', [0.00000000])])
    MIX[1] = dict([('name', 'fuel'),
                   ('SigmaT', 0.41666700),
                   ('Chi',    1.00000000),
                   ('SigmaS', 0.33400000),
                   ('NuSigF', 0.17800000)])
    forest.set_materials(ngroups=1, mix=MIX)
    forest_refl.set_materials(ngroups=1, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    #forest.memory()
    forest_refl.run()
    #plot_data(PROB_NAME+'.vtk')
    #forest.clean()

if __name__ == "__main__":
    test();
    #import doctest
    #doctest.testmod()

