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
    TESTDIR = "2dg7_c5g7/"
    PROB_NAME = "2dg7_c5g7_refl_logmesh"
    DESCRIPTION = ("3x3 pins representing a quarter.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 2])
    #METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 0), ('n_ref', 0)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', '2x2Pins2d'),
        ('nnodes', [3, 3]),
        ('length', [[21.42]*3, [21.42]*3]),
        ('components', [
            [4, 3, 0],
            [3, 4, 0],
            [1, 1, 2]]),
        ('boundary', [2, 0, 2, 0]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [10, 17]),
        ('components', [[0]*10 for i in range(17)]),
        ('length', [
            [0.3, 0.3, 0.6, 0.6, 1.2, 1.2, 1.6, 3.2, 5.4, 7.02], [1.26]*17,
        ]),
    ])
    LATT[1] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [17, 10]),
        ('components', [[0]*17 for i in range(10)]),
        ('length', [
            [1.26]*17, [0.3, 0.3, 0.6, 0.6, 1.2, 1.2, 1.6, 3.2, 5.4, 7.02],
        ]),
    ])
    LATT[2] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [10, 10]),
        ('components', [[0] * 10 for i in range(10)]),
        ('length', [
            [0.3, 0.3, 0.6, 0.6, 1.2, 1.2, 1.6, 3.2, 5.4, 7.02],
            [0.3, 0.3, 0.6, 0.6, 1.2, 1.2, 1.6, 3.2, 5.4, 7.02],
        ]),
    ])
    LATT[3] = dict([
        ('type', 'pin_map'),
        ('name', 'MOX'),
        ('nnodes', [17, 17]),
        ('components', [
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
            [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
            [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
            [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
            [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
            [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 6, 4, 4, 6, 4, 4, 5, 4, 4, 6, 4, 4, 6, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
            [2, 3, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 4, 4, 6, 3, 2],
            [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
            [2, 3, 3, 6, 3, 4, 4, 4, 4, 4, 4, 4, 3, 6, 3, 3, 2],
            [2, 3, 3, 3, 3, 6, 3, 3, 6, 3, 3, 6, 3, 3, 3, 3, 2],
            [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
            [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
        ]),
        ('water_gap', [0., 0., 0., 0.]),
        ('length', [[1.26]*17 for i in range(2)]),
    ])
    LATT[4] = dict([
        ('type', 'pin_map'),
        ('name', 'UOX'),
        ('nnodes', [17, 17]),
        ('components', [
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
            [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 6, 1, 1, 6, 1, 1, 5, 1, 1, 6, 1, 1, 6, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 6, 1, 1, 1],
            [1, 1, 1, 1, 1, 6, 1, 1, 6, 1, 1, 6, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        ]),
        ('water_gap', [0., 0., 0., 0.]),
        ('length', [[1.26]*17 for i in range(2)]),
    ])
    PIN = dict()
    PIN[0] = newpin(ptype='box', mat=[0], name='moderator');
    PIN[1] = newpin(ptype='pin', mat=[0, 1], fradius=0.54, name='UO2');
    PIN[2] = newpin(ptype='pin', mat=[0, 2], fradius=0.54, name='MOX4.3');
    PIN[3] = newpin(ptype='pin', mat=[0, 3], fradius=0.54, name='MOX7.0');
    PIN[4] = newpin(ptype='pin', mat=[0, 4], fradius=0.54, name='MOX8.7');
    PIN[5] = newpin(ptype='pin', mat=[0, 5], fradius=0.54, name='fisschamber');
    PIN[6] = newpin(ptype='pin', mat=[0, 6], fradius=0.54, name='guidetube');
    forest.set_geometry(dim=2, core=CORE, lattices=LATT, pins=PIN)
    # Materials ---------------------------------------------------------------
    from get_c5g7_materials import get_c5g7_materials
    MIX = get_c5g7_materials()
    forest.set_materials(ngroups=7, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    forest.clean()


