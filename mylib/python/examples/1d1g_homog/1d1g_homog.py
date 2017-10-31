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

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "1d1g_homog/"
    PROB_NAME = "1d1g_homog"
    DESCRIPTION = (
        'Homogeneous slab 10.0 cm (it can be a 2.0 cm or 10.0 cm slab) \n'
        'The reference for sn = 96 is keff = 0.662951 (2cm slab) \n'
        'The reference for sn = 96 is keff = 2.00193 (10cm slab) \n'
        'This problem is provided in the following paper: \n \n'
        'M. Capilla,  CF. Talavera,  D. Ginestar,  G. Verdu,  \n'
        '"A nodal collocation method for the calculation of \n'
        'the lambda modes of the P_L equations"; \n'
        'Annals of Nuclear Energy, 32,  (2005), 1825--1853.')
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport', dsa=False, quad=['GaussLegendre', 96])
    #METHOD = setmethod(mtype='transport', quad=['GaussLegendre', 96])
    FE_SETTINGS= dict([('degree', 4), ('n_ref', 2)])
    #ALGEBRA = setalgebra(matfree=True)
    ALGEBRA = setalgebra(matfree=False)
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    # Geometry ----------------------------------------------------------------
    CORE = dict([('composed', True),
                 ('name', 'mini core'),
                 ('nnodes', [1]),
                 ('length', [[2.0]]),
                 ('components', [0]),
                 ('boundary', [0, 0])])
    LATT = dict()
    LATT[0] = dict([('type', 'grid'),
                    ('name', 'assembly'),
                    ('nnodes', [5]),
                    ('components', [0, 0, 0, 0, 0]),
                    ('length', [0.4, 0.4, 0.4, 0.4, 0.4])])
    forest.set_geometry(dim=1, core=CORE, lattices=LATT)
    # Materials ---------------------------------------------------------------
    MIX = dict()
    MIX[0] = dict([('name', 'fuel'),
                   ('SigmaT', 1.0000000000),
                   ('Chi',    1.0000000000),
                   ('SigmaS', 0.9000000000),
                   ('NuSigF', 0.2500000000)])
    forest.set_materials(ngroups=1, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    forest.clean()


