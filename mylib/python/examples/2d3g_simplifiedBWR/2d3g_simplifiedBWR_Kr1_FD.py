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
from prob.utils_ploting_mat import plot_materials, plot_fluxes

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "2d3g_simplifiedBWR/"
    PROB_NAME = "2d3g_simplifiedBWR_Kr1_FD"
    DESCRIPTION = ("Lattice problem.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    ALGEBRA = setalgebra(matfree=True, eig=["Krylov:1:3",1.e-7,1000,"FD"], inner=["Krylov",1.e-9,1000,2])
    METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 4])
#    METHOD = setmethod(mtype='diffusion')

    FE_SETTINGS= dict([('degree', 2), ('n_ref', 1)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', 'Lattice'),
        ('nnodes', [1, 1]),
        ('length', [[14.0], [14.0]]),
        ('components', [0]),
        ('boundary', [0, 0, 0, 0]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [19, 19]),
        ('components', [
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]]),
        ('length', [
            [0.825, 1.15, 0.25, 1.15, 0.25, 1.15, 0.25, 1.15, 0.25, 1.15,
             0.25, 1.15, 0.25, 1.15, 0.25, 1.15, 0.25, 1.15, 0.825],
            [0.825, 1.15, 0.25, 1.15, 0.25, 1.15, 0.25, 1.15, 0.25, 1.15,
             0.25, 1.15, 0.25, 1.15, 0.25, 1.15, 0.25, 1.15, 0.825]]),
    ])
    forest.set_geometry(dim=2, core=CORE, lattices=LATT)
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
    #plot_materials(PROB_NAME)
    #plot_fluxes(PROB_NAME)
    #forest.clean()


