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
from prob.runforest import newpin
from prob.settings import setmethod
# we use numpy for the data
#import numpy as np

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "2d3g_lattice/"
    PROB_NAME = "2d3g_lattice"
    DESCRIPTION = ("Lattice problem.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    #METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 4])
    METHOD = setmethod(mtype='diffusion')
    forest.set_settings(method=METHOD)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', 'Lattice'),
        ('nnodes', [1, 1]),
        ('length', [[17.0], [17.0]]),
        ('components', [0]),
        ('boundary', [2, 2, 2, 2]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'pin_map'),
        ('name', 'assembly'),
        ('nnodes', [17, 17]),
        ('components', [[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]),
        ('water_gap', [0., 0., 0., 0.]),
        ('length', [[1.0]*17 for i in range(2)]),
    ])
    PIN = dict()
    PIN[0] = newpin(ptype='pin', mat=[0, 1], fradius=0.4, name='Pin');
    forest.set_geometry(dim=2, core=CORE, lattices=LATT, pins=PIN)
    # Materials ---------------------------------------------------------------
    MIX = dict()
    MIX[0] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.3333694980, 0.5887110656, 1.6512246291]),
        ('Chi',    [0.9996892490, 0.0003391680, 0.0000000000]),
        ('SigmaS', [[0.2432946408, 0.0000000000, 0.0000000000],
                    [0.0898364840, 0.4413753398, 0.0000203109],
                    [0.0000387911, 0.1465683257, 1.6300848232]]),
        ('NuSigF', [0.0000000000, 0.0000000000, 0.0000000000]),
    ])
    MIX[1] = dict([
        ('name', 'fuel'),
        ('SigmaT', [0.2822058997, 0.4997685502, 0.4323754911]),
        ('Chi',    [0.9996892490, 0.0003391680, 0.0000000000]),
        ('SigmaS', [[0.2760152893, 0.0000000000, 0.0000000000],
                    [0.0011230014, 0.4533430274, 0.0000378305],
                    [0.0000000000, 0.0014582502, 0.2823864370]]),
        ('NuSigF', [0.0076796277, 0.0234281944, 0.2734408326]),
    ])
    forest.set_materials(ngroups=3, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    forest.clean()


