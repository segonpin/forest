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
from prob.settings import setmethod, setalgebra
#from prob.utils import plot_materials, plot_fluxes
# we use numpy for the data
#import numpy as np

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "2d2g_c4/"
    PROB_NAME = "2d2g_c4"
    DESCRIPTION = ("C4 problem, obtained from: \n"
    "Cavarec, C., Perron, J.F., Verwaerde, D. and West, J.P., 1994.\n"
    "Benchmark calculations of power distribution within assemblies \n"
    "(No. NEA-NSC-DOC 94-28). Nuclear Energy Agency.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    #METHOD = setmethod(mtype='transport', dsa=False, quad=['LevelSymType2', 2])
    #METHOD = setmethod(mtype='transport', dsa=True, quad=['LevelSymType2', 2])
    METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 2])
    #METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 5), ('n_ref', 0)])
    #ALGEBRA = setalgebra(matfree=False)
    ALGEBRA = setalgebra(matfree=True)
    #OUTPUT = dict([("homogenize", True)])
    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', 'Lattice'),
        ('nnodes', [2, 2]),
        ('length', [[21.42, 21.42], [21.42, 21.42]]),
        ('components', [[0, 2], [2, 0]]),
        ('boundary', [2, 0, 2, 0]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [17, 17]),
        ('components', [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1],
                        [1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 5, 1, 1, 5, 1, 1, 7, 1, 1, 5, 1, 1, 5, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1],
                        [1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]),
        ('length', [[1.26]*17, [1.26]*17]),
    ])
    LATT[1] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [17, 17]),
        ('components', [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1],
                        [1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 8, 1, 1, 8, 1, 1, 7, 1, 1, 8, 1, 1, 8, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1],
                        [1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]),
        ('length', [[1.26]*17, [1.26]*17]),
    ])
    LATT[2] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [17, 17]),
        ('components', [[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2],
                        [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
                        [2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2],
                        [2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2],
                        [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
                        [2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2],
                        [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                        [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                        [2, 3, 5, 4, 4, 5, 4, 4, 7, 4, 4, 5, 4, 4, 5, 3, 2],
                        [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                        [2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2],
                        [2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2],
                        [2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2],
                        [2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2],
                        [2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2],
                        [2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2],
                        [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]]),
        ('length', [[1.26]*17, [1.26]*17]),
    ])
    LATT[3] = dict([
        ('type', 'grid'),
        ('name', 'assembly'),
        ('nnodes', [17, 17]),
        ('components', [[6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6],
                        [6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6]]),
        ('length', [[1.26]*17, [1.26]*17]),
    ])
    forest.set_geometry(dim=2, core=CORE, lattices=LATT)
    # Materials ---------------------------------------------------------------
    MIX = dict()
    MIX[1] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 0.833333333333333]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.247777777777778, 0.0000000000],
                    [0.020, 0.7333333333333333]]),
        ('NuSigF', [0.0050, 0.125]),
    ])
    MIX[2] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 0.833333333333333]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.247777777777778, 0.0000000000],
                    [0.015, 0.6333333333333333]]),
        ('NuSigF', [0.0075, 0.300]),
    ])
    MIX[3] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 0.833333333333333]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.247777777777778, 0.0000000000],
                    [0.015, 0.5833333333333333]]),
        ('NuSigF', [0.0075, 0.375]),
    ])
    MIX[4] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 0.833333333333333]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.247777777777778, 0.0000000000],
                    [0.015, 0.5333333333333332]]),
        ('NuSigF', [0.0075, 0.450]),
    ])
    MIX[5] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 0.833333333333333]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.251777777777778, 0.0000000000],
                    [0.025, 0.8133333333333332]]),
        ('NuSigF', [0.000, 0.000]),
    ])
    MIX[6] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 1.666666666666667]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.226777777777778, 0.0000000000],
                    [0.050, 1.6266666666666665]]),
        ('NuSigF', [0.000, 0.000]),
    ])
    MIX[7] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 0.833333333333333]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.251777777777778, 0.0000000000],
                    [0.025, 0.8133333333333332]]),
        ('NuSigF', [1e-7, 3e-6]),
    ])
    MIX[8] = dict([
        ('name', 'moderator'),
        ('SigmaT', [0.277777777777778, 0.833333333333333]),
        ('Chi',    [1.000000000000000, 0.000000000000000]),
        ('SigmaS', [[0.227777777777778, 0.0000000000],
                    [0.010, 0.0333333333333332]]),
        ('NuSigF', [0.000, 0.000]),
    ])
    forest.set_materials(ngroups=2, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    #forest.clean()
    #plot_materials(PROB_NAME)
    #plot_fluxes(PROB_NAME)


