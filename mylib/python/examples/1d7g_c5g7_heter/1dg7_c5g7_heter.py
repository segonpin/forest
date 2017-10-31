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
from prob.runforest import RunForest, newpin
from prob.settings import setmethod, setalgebra

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "1dg7_c5g7_heter/"
    PROB_NAME = "1dg7_c5g7_heter"
    DESCRIPTION = ("One-dimensional heterogeneous version of the c5g7.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport', quad=['GaussLegendre', 48])
    #METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 0), ('n_ref', 0)])
    #ALGEBRA = setalgebra(matfree=True, eig=["PI",1.e-7,1000,"Standard"], inner=["Krylov",1.e-9,1000])
    #ALGEBRA = setalgebra(matfree=True, eig=["Krylov:2",1.e-7,1000,"Standard"], inner=["Krylov",1.e-9,1000, 7])
    #ALGEBRA = setalgebra(matfree=True, eig=["PI",1.e-7,1000, "Standard"], inner=["Krylov",1.e-9,1000,7])
    ALGEBRA = setalgebra(matfree=True, form = "FD", eig=["PI",1.e-7,1000], mg=["GS",1.e-7,10], inner=["Krylov",1.e-9,1000])
    #ALGEBRA = setalgebra(matfree=True, form = "FD", eig=["PI",1.e-7,1000], mg=["krylov",1.e-7,100], inner=["Krylov",1.e-9,1000])
    #ALGEBRA = setalgebra(matfree=True, eig=["Krylov:1:3",1.e-7,1000, "Standard"], inner=["Krylov",1.e-9,1000,10])
    #ALGEBRA = setalgebra(matfree=True, eig=["Krylov:1:3",1.e-7,1000, "FD"], inner=["Krylov",1.e-9,1000,10])

    forest.set_settings(method=METHOD, fe_settings=FE_SETTINGS, algebra=ALGEBRA)
    # Geometry ----------------------------------------------------------------
    CORE = dict([
        ('composed', True),
        ('name', 'Core'),
        ('nnodes', [3]),
        ('length', [[21.42, 21.42, 21.42]]),
        ('components', [[2, 1, 2]]),
        ('boundary', [2, 2]),
    ])
    LATT = dict()
    LATT[1] = dict([
        ('type', 'pin_map'),
        ('name', 'MOX'),
        ('nnodes', [17]),
        ('components', [[2, 3, 6, 4, 4, 6, 4, 4, 5, 4, 4, 6, 4, 4, 6, 3, 2]]),
        #('components', [[2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2]]),
        ('water_gap', [0., 0.]),
        ('length', [[1.26]*17]),
    ])
    LATT[2] = dict([
        ('type', 'pin_map'),
        ('name', 'UOX'),
        ('nnodes', [17]),
        ('components', [[1, 1, 6, 1, 1, 6, 1, 1, 5, 1, 1, 6, 1, 1, 6, 1, 1]]),
        #('components', [[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]),
        ('water_gap', [0., 0.]),
        ('length', [[1.26]*17]),
    ])
    PIN = dict()
    PIN[0] = newpin(ptype='box', mat=[0], name='moderator');
    PIN[1] = newpin(ptype='pin', mat=[0, 1], fradius=0.54, name='UO2');
    PIN[2] = newpin(ptype='pin', mat=[0, 2], fradius=0.54, name='MOX4.3');
    PIN[3] = newpin(ptype='pin', mat=[0, 3], fradius=0.54, name='MOX7.0');
    PIN[4] = newpin(ptype='pin', mat=[0, 4], fradius=0.54, name='MOX8.7');
    PIN[5] = newpin(ptype='pin', mat=[0, 5], fradius=0.54, name='fisschamber');
    PIN[6] = newpin(ptype='pin', mat=[0, 6], fradius=0.54, name='guidetube');
    forest.set_geometry(dim=1, core=CORE, lattices=LATT, pins=PIN)
    # Materials ---------------------------------------------------------------
    from get_c5g7_materials import get_c5g7_materials
    MIX = get_c5g7_materials()
    forest.set_materials(ngroups=7, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    #forest.clean()


