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
from prob.runforest import RunForest
from prob.settings import setmethod, setalgebra

#import numpy as np
#

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "2d1g_heter/"
    PROB_NAME = "2d1g_heter"
    DESCRIPTION = (
        "THeterogeneous square of 7x7 regions. \n"
        "The reference results have been obtained with dragon, usging \n"
        "degree three (d:=3), no refinement (s:=1) and type 2 \n"
        "quadrature set (QUAD = 2); \n"
        "The reference for sn = 02 is keff = 5.755775E-01 \n"
        "The reference for sn = 04 is keff = 5.806699E-01 \n"
        "The reference for sn = 06 is keff = 5.900330E-01 \n"
        "The reference for sn = 08 is keff = 5.941040E-01 \n"
        "The reference for sn = 12 is keff = 5.957445E-01 \n"
        "The reference for sn = 16 is keff = 5.963024E-01 \n"
        "The reference for sn = 20 is keff = 5.964814E-01 \n \n"
        "This problem is a generalization of the 1D seven \n"
        "region problem provided by: \n \n"
        "M. Capilla,  CF. Talavera,  D. Ginestar,  G. Verdu, \n"
        "'A nodal collocation method for the calculation of \n"
        "the lambda modes of the P_L equations'; \n"
        "Annals of Nuclear Energy, 32,  (2005), 1825-1853.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport',dsa=False, quad=['LevelSymType2', 4])
    #METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 2])
    #METHOD = setmethod(mtype='diffusion')
    ALGEBRA = setalgebra(matfree=True)
    FE_SETTINGS= dict([('degree', 3), ('n_ref', 1)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS,algebra=ALGEBRA)
    # Geometry ----------------------------------------------------------------
    CORE = dict([('composed', True),
                 ('name', 'mini core'),
                 ('nnodes', [1, 1]),
                 ('length', [[18.0], [18.0]]),
                 ('components', [0]),
                 ('boundary', [0, 0, 0, 0])])
    LATT = dict()
    LATT[0] = dict([('type', 'grid'),
                    ('name', 'assembly'),
                    ('nnodes', [8, 8]),
                    ('components', [[0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 1, 0, 1, 1, 0, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 1, 0, 1, 1, 0, 1, 0],
                                    [0, 1, 0, 1, 1, 0, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0],
                                    [0, 1, 0, 1, 1, 0, 1, 0],
                                    [0, 0, 0, 0, 0, 0, 0, 0]]),
                    ('length', [[2.7, 2.4, 2.7, 1.2, 1.2, 2.7, 2.4, 2.7],
                                [2.7, 2.4, 2.7, 1.2, 1.2, 2.7, 2.4, 2.7]])])
    forest.set_geometry(dim=2, core=CORE, lattices=LATT)
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
    # Running forest ----------------------------------------------------------
    forest.run()
    #forest.clean()



