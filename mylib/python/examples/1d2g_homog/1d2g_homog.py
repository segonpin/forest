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
import numpy as np
#
from prob.runforest import RunForest
from prob.settings import setmethod

if __name__ == "__main__":

    # -------------------------------------------------------------------------
    TESTDIR = "1d2g_homog/"
    PROB_NAME = "1d2g_homog"
    DESCRIPTION = (
        "Here we define the different materials associated to the \n"
        "material id specified when defining the geometry  \n"
        "The diffusion cross sections are: \n"
        "#   D1     D2    SIGMA_A1  SIGMA_A2 SIGMA_S12 SIGMA_F1   SIGMA_F2  \n"
        "1.51244 0.297925 0.00846854 0.0638224 0.0190582 0.00620569 0.10209, \n"
        "so then the transport cross sections should be   \n"
        "# SIGMA_T1    SIGMA_T2     SIGMA_S1     SIGMA_S2     SIGMA_S12 SIGMA_F1   SIGMA_F2 \n"
        "0.2203944178  1.1188498224 0.1928676778 1.0550274224 0.0190582 0.00620569 0.10209, \n"
        "And the reference eigenvalue should be lambda = 1.332924e+00 \n"
        "sa1 =  0.0084685 \n"
        "sa2 =  0.063822 \n"
        "s12 =  0.019058 \n"
        "sf1 =  0.0062057 \n"
        "sf2 =  0.10209 \n"
        "octave:5> (sf1*sa2 + sf2*s12)/((sa1+s12)*sa2+sa2) \n"
        "ans =  1.33292550586827 ?? \n"
        "#the formulas for the cross sections are \n"
        "SigmaT = 1/(3*D) \n"
        "SigmaS_{g,g} = SigmaT - SigmaA_{g->g} + Sum_{h not g}(SigmaS_{g->h}) \n"
        "I am not sure now, but it depends on the scattering in and out, as \n"
        "they explain here: \n"
        "Reactor Physics: Multigroup Diffusion  \n"
        "prepared by W m . J. Garland, Professor, \n"
        "Department of Engineering Physics,   \n"
        "McMaster University, Hamilton, Ontario, Canada \n"
        "http://www.nuceng.ca/ep4d3/text/9-multigroup.pdf.")
    forest = RunForest(prob_name=PROB_NAME, description=DESCRIPTION)
    # Settings ----------------------------------------------------------------
    METHOD = setmethod(mtype='diffusion')
    forest.set_settings(method=METHOD)
    # Geometry ----------------------------------------------------------------
    CORE = dict([('composed', True),
                 ('name', 'GeorgiaTech1D'),
                 ('nnodes', [1]),
                 ('length', [[15.24]]),
                 ('components', [0]),
                 ('boundary', [2, 2])])
    LATT = dict()
    LATT[0] = dict([('type', 'grid'),
                    ('name', 'assembly'),
                    ('nnodes', [6]),
                    ('components', [0, 0, 0, 0, 0, 0]),
                    ('length', [2.54, 2.54, 2.54, 2.54, 2.54, 2.54])])
    forest.set_geometry(dim=1, core=CORE, lattices=LATT)
    # Materials ---------------------------------------------------------------
    MIX = dict()
    MIX[0] = dict([('name', 'moderator'),
                   ('SigmaT', np.array([0.2203944178, 1.1188498224])),
                   ('Chi',    np.array([1.0000000000, 0.0000000000])),
                   ('SigmaS', np.array([[0.1928676778, 0.0000000000],
                                        [0.0190582000, 1.0550274224]])),
                   ('NuSigF', np.array([0.0062056900, 0.1020900000]))])
    forest.set_materials(ngroups=2, mix=MIX)
    # Running forest ----------------------------------------------------------
    forest.run()
    forest.clean()

