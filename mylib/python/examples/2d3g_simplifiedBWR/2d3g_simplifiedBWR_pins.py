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
from prob.settings import setmethod, setalgebra

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "2d3g_simplifiedBWR/"
    PROB_NAME = "2d3g_simplifiedBWR_pins"
    DESCRIPTION = ("Lattice problem.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 2])
    #METHOD = setmethod(mtype='diffusion')
    ALGEBRA = setalgebra(matfree=True, eig=["Krylov:3:7",1.e-7,1000], inner=["Krylov",1.e-8,1000,2])
    #ALGEBRA = setalgebra(matfree=True, eig=["Krylov:1",1.e-7,1000])
    #ALGEBRA = setalgebra(matfree=True, eig=["PI",1.e-7,1000])

    forest.set_settings(method=METHOD, algebra=ALGEBRA)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', 'Lattice'),
        ('nnodes', [3, 3]),
        ('length', [[0.7, 12.6, 0.7], [0.7, 12.6, 0.7]]),
        ('components', [
            [3, 1, 3],
            [2, 0, 2],
            [3, 1, 3]]),
        ('boundary', [2, 2, 2, 2]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'pin_map'),
        ('name', 'assembly'),
        ('nnodes', [9, 9]),
        ('components', [
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 0, 0, 1, 1, 1],
            [1, 1, 1, 0, 0, 0, 1, 1, 1],
            [1, 1, 1, 0, 0, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1],
            [1, 1, 1, 1, 1, 1, 1, 1, 1]]),
        ('water_gap', [0.0, 0.0, 0.0, 0.0]),
        ('length', [[1.4]*9 for i in range(2)]),
    ])
    LATT[1] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [9, 1]),
        ('components', [0, 0, 0, 0, 0, 0, 0, 0, 0]),
        ('length', [[1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4], [0.7]]),
    ])
    LATT[2] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [1, 9]),
        ('components', [[0], [0], [0], [0], [0], [0], [0], [0], [0]]),
        ('length', [[0.7], [1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4]]),
    ])
    LATT[3] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [1, 1]),
        ('components', [0]),
        ('length', [[0.7], [0.7]]),
    ])
    PIN = dict()
    PIN[0] = newpin(ptype='box', mat=[0], name='Pin');
    PIN[1] = newpin(ptype='pinbox', mat=[0, 1], fradius=0.575, name='Pin');
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


