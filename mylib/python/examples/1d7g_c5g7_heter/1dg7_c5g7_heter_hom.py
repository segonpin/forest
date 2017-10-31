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
# extra functions to make our life easier
from prob.runforest import newpin
from prob.utils import combine_xml, combine_names
from prob.settings import setmethod

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    #TESTDIR = "1dg7_c5g7_heter/"
    #PROB_NAME = "1dg7_c5g7_heter"
    DESCRIPTION = ("One-dimensional heterogeneous version of the c5g7.")
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='transport', quad=['GaussLegendre', 2])
    #METHOD = setmethod(mtype='diffusion')
    FE_SETTINGS= dict([('degree', 0), ('n_ref', 0)])
    # Materials ---------------------------------------------------------------
    from get_c5g7_materials import get_c5g7_materials
    MIX = get_c5g7_materials()
    # Geometry ----------------------------------------------------------------
    PIN = dict()
    PIN[0] = newpin(ptype='box', mat=[0], name='moderator');
    PIN[1] = newpin(ptype='pin', mat=[0, 1], fradius=0.54, name='UO2');
    PIN[2] = newpin(ptype='pin', mat=[0, 2], fradius=0.54, name='MOX4.3');
    PIN[3] = newpin(ptype='pin', mat=[0, 3], fradius=0.54, name='MOX7.0');
    PIN[4] = newpin(ptype='pin', mat=[0, 4], fradius=0.54, name='MOX8.7');
    PIN[5] = newpin(ptype='pin', mat=[0, 5], fradius=0.54, name='fisschamber');
    PIN[6] = newpin(ptype='pin', mat=[0, 6], fradius=0.54, name='guidetube');
    LATT = dict()
    LATT[0] = dict([
        ('type', 'grid'),
        ('name', 'REFLECTOR'),
        ('nnodes', [17]),
        ('components', [0]*17 ),
        ('length', [[1.26]*17]),
    ])
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
    COREASS1 = dict([
        ('composed', True),
        ('name', 'Core'),
        ('nnodes', [1]),
        ('length', [[21.42]]),
        ('components', [[1]]),
        ('boundary', [2, 2]),
    ])
    COREASS2 = dict([
        ('composed', True),
        ('name', 'Core'),
        ('nnodes', [1]),
        ('length', [[21.42]]),
        ('components', [[2]]),
        ('boundary', [2, 2]),
    ])
    # Renaming files ----------------------------------------------------------
    PROB_ASS1 = '1dg7_c5g7_heter_ASS1'
    PROB_ASS2 = '1dg7_c5g7_heter_ASS2'
    # -------------------------------------------------------------------------
    probAss1 = RunForest(prob_name=PROB_ASS1, description=DESCRIPTION)
    probAss1.set_materials(ngroups=7, mix=MIX)
    probAss1.set_geometry(dim=1, core=COREASS1, lattices=LATT, pins=PIN)
    OUTPUT = dict([("homogenize", True)])
    probAss1.set_settings(method=METHOD, fe_settings=FE_SETTINGS, output=OUTPUT)
    # -------------------------------------------------------------------------
    probAss2 = RunForest(prob_name=PROB_ASS2, description=DESCRIPTION)
    probAss2.set_materials(ngroups=7, mix=MIX)
    probAss2.set_geometry(dim=1, core=COREASS2, lattices=LATT, pins=PIN)
    probAss2.set_settings(method=METHOD, fe_settings=FE_SETTINGS, output=OUTPUT)
    # Running forest ----------------------------------------------------------
    # -------------------------------------------------------------------------
    probAss1.run()
    probAss1.clean()
    # -------------------------------------------------------------------------
    probAss2.run()
    probAss2.clean()

    # -------------------------------------------------------------------------
    # Assembly homogenized
    # -------------------------------------------------------------------------
    PROB_HOMASS = '1dg7_c5g7_heter_HOMASS'
    DESCRIPTION = ("Homogenized version of the c5g7.")
    probHomAss = RunForest(prob_name=PROB_HOMASS, description=DESCRIPTION)
    # Geometry ----------------------------------------------------------------
    COREHOM = dict([
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
        ('components', [[0]*17]),
        ('water_gap', [0., 0.]),
        ('length', [[1.26]*17]),
    ])
    LATT[2] = dict([
        ('type', 'pin_map'),
        ('name', 'UOX'),
        ('nnodes', [17]),
        ('components', [[1]*17]),
        ('length', [[1.26]*17]),
    ])
    PIN = dict()
    for i in range(2):
        PIN[i] = newpin(ptype='pin', mat=[i, i], fradius=0.54, name='mix'+str(i));
    probHomAss.set_geometry(dim=1, core=COREHOM, lattices=LATT, pins=PIN)
    # Materials ---------------------------------------------------------------
    ass_list = [PROB_ASS1,PROB_ASS2]
    ass_list = [ass+".ass.mat.xml" for ass in ass_list]

    #print combine_names(ass_list,".ass.mat.xml")
    #combine_xml(ass_list).write(combine_names(ass_list,".ass.mat.xml"))

    probHomAss.attach_materials(combine_xml(ass_list))
    # Settings ----------------------------------------------------------------
    probHomAss.set_settings(method=METHOD, fe_settings=FE_SETTINGS)
    ## Running forest ---------------------------------------------------------
    probHomAss.run()
    probHomAss.clean()

    # -------------------------------------------------------------------------
    # Pin homogenized
    # -------------------------------------------------------------------------
    PROB_HOMPIN = '1dg7_c5g7_heter_HOMPIN'
    DESCRIPTION = ("Homogenized version of the c5g7.")
    probHomPin = RunForest(prob_name=PROB_HOMPIN, description=DESCRIPTION)
    # Geometry ----------------------------------------------------------------
    COREHOM = dict([
        ('composed', True),
        ('name', 'Core'),
        ('nnodes', [3]),
        ('length', [[21.42, 21.42, 21.42]]),
        ('components', [[0, 1, 2]]),
        ('boundary', [2, 2]),
    ])
    LATT = dict()
    LATT[0] = dict([
        ('type', 'pin_map'),
        ('name', 'UOX'),
        ('nnodes', [17]),
        ('components', range(0,17) ),
        ('length', [[1.26]*17]),
    ])
    LATT[1] = dict([
        ('type', 'pin_map'),
        ('name', 'MOX'),
        ('nnodes', [17]),
        ('components', [range(17,34)]),
        ('length', [[1.26]*17]),
    ])
    LATT[2] = dict([
        ('type', 'pin_map'),
        ('name', 'UOX'),
        ('nnodes', [17]),
        ('components', [range(34,51)]),
        ('length', [[1.26]*17]),
    ])
    PIN = dict()
    for i in range(51):
        PIN[i] = newpin(ptype='pin', mat=[i, i], fradius=0.54, name='mix'+str(i));
    probHomPin.set_geometry(dim=1, core=COREHOM, lattices=LATT, pins=PIN)
    # Materials ---------------------------------------------------------------
    ass_list = [PROB_ASS2,PROB_ASS1,PROB_ASS2]
    ass_list = [ass+".pin.mat.xml" for ass in ass_list]
    probHomPin.attach_materials(combine_xml(ass_list))
    # Settings ----------------------------------------------------------------
    probHomPin.set_settings(method=METHOD, fe_settings=FE_SETTINGS)
    ## Running forest ----------------------------------------------------------
    probHomPin.run()
    probHomPin.clean()

