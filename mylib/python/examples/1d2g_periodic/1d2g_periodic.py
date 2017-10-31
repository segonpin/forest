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
    TESTDIR = "1d2g_periodic"
    PROB_HOM = "1d2g_periodic_homog"
    PROB_HET = "1d2g_periodic_heter"
    DESCRIPTION = ("Simplified version of the georgiatech.")
    forest_hom = RunForest(prob_name=PROB_HOM, description=DESCRIPTION)
    forest_het = RunForest(prob_name=PROB_HET, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport', quad=['GaussLegendre', 48])
    #METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 2), ('n_ref', 2)])
    OUTPUT = dict([("homogenize", False)])
    forest_hom.set_settings(method=METHOD, fe_settings=FE_SETTINGS, output=OUTPUT)
    forest_het.set_settings(method=METHOD, fe_settings=FE_SETTINGS, output=OUTPUT)
    # Geometry ----------------------------------------------------------------
    CORE_HOM = dict([
        ('composed', True),
        ('name', 'Core'),
        ('nnodes', [2]),
        ('length', [[15.24, 15.24]]),
        ('components', [[0, 0]]),
        ('boundary', [2, 2]),
    ])
    CORE_HET = dict([
        ('composed', True),
        ('name', 'Core'),
        ('nnodes', [2]),
        ('length', [[15.24, 15.24]]),
        ('components', [[1, 1]]),
        ('boundary', [2, 2]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'grid'),
        ('name', 'homog ass'),
        ('nnodes', [6]),
        ('components', [3]*6 ),
        ('length', [[2.54]*6]),
    ])
    LATT[1] = dict([
        ('type', 'grid'),
        ('name', 'heter ass'),
        ('nnodes', [6]),
        ('components', [[0, 3, 3, 3, 3, 0]]),
        ('length', [1.158, 3.231, 3.231, 3.231, 3.231, 1.158]),
    ])
    forest_hom.set_geometry(dim=1, core=CORE_HOM, lattices=LATT)
    forest_het.set_geometry(dim=1, core=CORE_HET, lattices=LATT)
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
    forest_hom.set_materials(ngroups=2, mix=MIX)
    forest_het.set_materials(ngroups=2, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest_hom.run()
    forest_hom.clean()

    forest_het.run()
    forest_het.clean()


