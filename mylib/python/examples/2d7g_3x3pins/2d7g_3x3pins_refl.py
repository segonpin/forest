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
    TESTDIR = "2dg7_3x3pins/"
    PROB_NAME = "2dg7_3x3pins_refl"
    DESCRIPTION = ("3x3 pins representing a quarter.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    #METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 4])
    METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 2), ('n_ref', 2)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', '2x2Pins2d'),
        ('nnodes', [1, 1]),
        ('length', [[3.78], [3.78]]),
        ('components', [0]),
        ('boundary', [2, 0, 2, 0]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'pin_map'),
        ('name', 'assembly'),
        ('nnodes', [3, 3]),
        ('components', [[1, 1, 0],
                        [2, 1, 0],
                        [0, 0, 0]]),
        ('length', [[1.26]*3 for i in range(2)]),
    ])
    PIN = dict()
    PIN[0] = newpin(ptype='box', mat=[0], name='moderator');
    PIN[1] = newpin(ptype='pin', mat=[0, 2], fradius=0.54, name='Pin');
    PIN[2] = newpin(ptype='pin', mat=[0, 1], fradius=0.34, name='Pin');
    forest.set_geometry(dim=2, core=CORE, lattices=LATT, pins=PIN)
    # Materials ---------------------------------------------------------------
    from get_c5g7_materials import get_c5g7_materials
    MIX = get_c5g7_materials()
    forest.set_materials(ngroups=7, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    forest.clean()


