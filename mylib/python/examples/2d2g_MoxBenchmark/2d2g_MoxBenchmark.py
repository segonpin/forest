# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:22:51 2016

@author: segonpin
"""
# we add the root directory for the python library
# if we do not have that directory in our path
import sys
sys.path.insert(0, '../..')
#
from prob.runforest import RunForest
from prob.settings import setmethod, setalgebra
#from prob.utils import plot_materials, plot_fluxes
# we use numpy for the data
#import numpy as np

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "2d2g_MoxBenchmark/"
    PROB_NAME = "2d2g_MoxBenchmark"
    DESCRIPTION = (
    'Mox Benchamark obtained from: \n'
    'M. Capilla, C. F. Talavera, D. Ginestar, and G Verdu. \n'
    'A nodal collocation approximation for the multi-dimensional PL equations\n'
    '2D applications. Annals of Nuclear Energy, 35(10):1820-1830, 2008.')
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    #METHOD = setmethod(mtype='transport', dsa=False, quad=['LevelSymType2', 2])
    #METHOD = setmethod(mtype='transport', dsa=True, quad=['LevelSymType2', 2])
    METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 2])
    #METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 1), ('n_ref', 1)])
    #ALGEBRA = setalgebra(matfree=False)
    ALGEBRA = setalgebra(matfree=True, eig=["PI",1.e-8, 10000])
    #ALGEBRA = setalgebra(matfree=True, eig=["PI",1.e-8, 10000,"Standard"])
    #ALGEBRA = setalgebra(matfree=True, eig=["Krylov",1.e-8,10000,"FD"])
    #ALGEBRA = setalgebra(matfree=True, eig=["Krylov",1.e-8,10000,"Standard"])
    #OUTPUT = dict([("homogenize", True)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    
    # Geometry ----------------------------------------------------------------   
    CORE = dict([('composed', True),
             ('name', 'mini core'),
             ('nnodes', [1, 1]),
             ('length', [[192.78], [192.78]]),
             ('components', [0]),
             ('boundary', [0, 0, 0, 0])])
    
    LATT = dict()
    LATT[0] = dict([('type', 'grid'),
                    ('name', 'assembly'),
                    ('nnodes', [9, 9]),
                    ('components', [[2, 2, 2, 2, 2, 2, 2, 2, 2],
                                    [2, 1, 1, 1, 1, 1, 1, 1, 2],
                                    [2, 1, 1, 1, 0, 1, 1, 1, 2],
                                    [2, 1, 1, 0, 1, 0, 1, 1, 2],
                                    [2, 1, 0, 1, 0, 1, 0, 1, 2],
                                    [2, 1, 1, 0, 1, 0, 1, 1, 2],
                                    [2, 1, 1, 1, 0, 1, 1, 1, 2],
                                    [2, 1, 1, 1, 1, 1, 1, 1, 2],
                                    [2, 2, 2, 2, 2, 2, 2, 2, 2]]),
                    ('length', [[21.42]*9, [21.42]*9])])
       
    forest.set_geometry(dim=2, core=CORE, lattices=LATT)
    # Materials ---------------------------------------------------------------
    MIX = dict()
    MIX[0] = dict([
        ('name', 'MoxFuel'),
        ('SigmaT',  [0.550, 1.060]),
        ('Chi',     [1.000, 0.000]),
        ('SigmaS', [[0.520, 0.000],
                    [0.015, 0.760]]),
        ('NuSigF', [0.0075, 0.450])])
    MIX[1] = dict([
        ('name', 'UO2'),
        ('SigmaT',  [0.570, 1.100]),
        ('Chi',     [1.000, 0.000]),
        ('SigmaS', [[0.540, 0.000],
                    [0.020, 1.000]]),
        ('NuSigF',  [0.005, 0.125])])
    MIX[2] = dict([
        ('name', 'Reflector'),
        ('SigmaT',  [0.611, 2.340]),
        ('Chi',     [1.000, 0.000]),
        ('SigmaS', [[0.560, 0.000],
                    [0.050, 2.300]]),
        ('NuSigF',  [0.000, 0.000])])
    
    forest.set_materials(ngroups=2, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    #forest.clean()
    #plot_materials(PROB_NAME)
    #plot_fluxes(PROB_NAME)

