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
from prob.runforest import newpin
from prob.settings import setmethod

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "2dg7_1x1pins/"
    PROB_NAME = "2dg7_1x1pins"
    DESCRIPTION = ("Lattice problem.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 4])
    #METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 0), ('n_ref', 1)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', 'Pin2d'),
        ('nnodes', [1, 1]),
        ('length', [[1.26], [1.26]]),
        ('components', [0]),
        ('boundary', [0, 0, 0, 0]),
        #('boundary', [2, 2, 2, 2]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'pin_map'),
        ('name', 'assembly'),
        ('nnodes', [1, 1]),
        ('components', [[0]]),
        ('water_gap', [0., 0., 0., 0.]),
        ('length', [[1.26], [1.26]]),
    ])
    PIN = dict()
    PIN[0] = newpin(ptype='pinnew', mat=[0, 2], fradius=0.54, name='Pin');
    forest.set_geometry(dim=2, core=CORE, lattices=LATT, pins=PIN)
    # Materials ---------------------------------------------------------------
    from get_c5g7_materials import get_c5g7_materials
    MIX = get_c5g7_materials()
    forest.set_materials(ngroups=7, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    #forest.clean()


