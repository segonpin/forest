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
from prob.settings import setmethod

if __name__ == "__main__":
    # -------------------------------------------------------------------------
    TESTDIR = "1d2g_georgiatech/"
    PROB_NAME = "1d2g_georgiatech"
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
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport', quad=['GaussLegendre', 2])
    FE_SETTINGS= dict([('degree', 3), ('n_ref', 0)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS)
    # Geometry ----------------------------------------------------------------
    CORE = dict([('composed', True),
                 ('name', 'GeorgiaTech'),
                 ('nnodes', [7]),
                 ('length', [[15.24]*7]),
                 #('components', [1, 2, 1, 2, 1, 2, 1]), # conf. 1
                 #('components', [1, 3, 1, 3, 1, 3, 1]), # conf. 2
                 ('components', [1, 4, 1, 4, 1, 4, 1]), # conf. 3
                 ('boundary', [2, 2])])
    LATT = dict()
    LATT[0] = dict([('type', 'grid'),
                    ('name', 'Type0'),
                    ('nnodes', [6]),
                    ('components', [1]*6),
                    ('length', [2.54]*6)])
    LATT[1] = dict([('type', 'grid'),
                    ('name', 'Type1'),
                    ('nnodes', [6]),
                    ('components', [0, 1, 2, 2, 1, 0]),
                    ('length', [1.158, 3.231, 3.231, 3.231, 3.231, 1.158])])
    LATT[2] = dict([('type', 'grid'),
                    ('name', 'Type2'),
                    ('nnodes', [6]),
                    ('components', [0, 1, 1, 1, 1, 0]),
                    ('length', [1.158, 3.231, 3.231, 3.231, 3.231, 1.158])])
    LATT[3] = dict([('type', 'grid'),
                    ('name', 'Type3'),
                    ('nnodes', [6]),
                    ('components', [0, 1, 3, 3, 1, 0]),
                    ('length', [1.158, 3.231, 3.231, 3.231, 3.231, 1.158])])
    LATT[4] = dict([('type', 'grid'),
                    ('name', 'Type4'),
                    ('nnodes', [6]),
                    ('components', [0, 3, 3, 3, 3, 0]),
                    ('length', [1.158, 3.231, 3.231, 3.231, 3.231, 1.158])])
    forest.set_geometry(dim=1, core=CORE, lattices=LATT)
    # Materials ---------------------------------------------------------------
    MIX = dict()
    MIX[0] = dict([('name', 'water'),
                   ('SigmaT', [0.1889751875, 1.4632718759]),
                   ('Chi',    [1.0000000000, 0.0000000000]),
                   ('SigmaS', [[0.1506751875, 0.0000000000],
                               [0.0380000000, 1.4535718759]]),
                   ('NuSigF', [0.0000000000, 0.0000000000])])
    MIX[1] = dict([('name', 'FuelI'),
                   ('SigmaT', [0.2262955419, 1.0119409026]),
                   ('Chi',    [1.0000000000, 0.0000000000]),
                   ('SigmaS', [[0.2005955419, 0.0000000000],
                               [0.0161000000, 0.9355409026]]),
                   ('NuSigF', [0.0067000000, 0.1241000000])])
    MIX[2] = dict([('name', 'FuelII'),
                   ('SigmaT', [0.2251643699, 0.9914733293]),
                   ('Chi',    [1.0000000000, 0.0000000000]),
                   ('SigmaS', [[0.1994643699, 0.0000000000],
                               [0.0156000000, 0.9013733293]]),
                   ('NuSigF', [0.0078000000, 0.1542000000])])
    MIX[3] = dict([('name', 'FuelIIg'),
                   ('SigmaT', [0.2172685004, 1.0605578534]),
                   ('Chi',    [1.0000000000, 0.0000000000]),
                   ('SigmaS', [[0.1901685004, 0.0000000000],
                               [0.0136000000, 0.5732578534]]),
                   ('NuSigF', [0.0056000000, 0.0187000000])])
    forest.set_materials(ngroups=2, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    forest.clean()



