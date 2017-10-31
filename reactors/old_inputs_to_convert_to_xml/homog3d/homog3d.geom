# ##########################################################
# Homogeneous square
# The reference results have been obtained with dragon,
# using degree three (d:=3),  and no refinement (s:=1)
# for two different quadrature set:
# * type 1 quadrature set (QUAD = 1):
#   The reference for sn = 02 is keff = 1.654333E+00
#   The reference for sn = 04 is keff = 1.672745E+00
#   The reference for sn = 06 is keff = 1.674243E+00
#   The reference for sn = 08 is keff = 1.674910E+00
#   The reference for sn = 12 is keff = 1.675369E+00
#   The reference for sn = 16 is keff = 1.675521E+00
#
# * type 2 quadrature set (QUAD = 2):
#   The reference for sn = 02 is keff = 1.654333E+00
#   The reference for sn = 04 is keff = 1.672745E+00
#   The reference for sn = 06 is keff = 1.674243E+00
#   The reference for sn = 08 is keff = 1.674910E+00
#   The reference for sn = 12 is keff = 1.675369E+00
#   The reference for sn = 16 is keff = 1.675521E+00
#
# This problem is a generalization of the 1D seven
# region problem provided by:
#
# M. Capilla,  CF. Talavera,  D. Ginestar,  G. Verdu, 
# "A nodal collocation method for the calculation of
# the lambda modes of the P_L equations";
# Annals of Nuclear Energy, 32,  (2005), 1825--1853.
# ##########################################################


# Definitions of materials

n_nodes 5 5 2
n_xsecs 1 

Materials 2
axial level 0
0   0   0   0   0 
0   0   0   0   0 
0   0   0   0   0 
0   0   0   0   0 
0   0   0   0   0 
axial level 1
0   0   0   0   0 
0   0   0   0   0 
0   0   0   0   0 
0   0   0   0   0 
0   0   0   0   0 


# I have to change that to have
# step_sizes
# x_lenght 2.7 2.4 2.7 2.4 2.7
# y_lenght 2.7 2.4 2.7 2.4 2.7
# z_lenght 1.0 1.0 1.0 1.0 1.0

step_length 
2.0 2.0 2.0 2.0 2.0 
2.0 2.0 2.0 2.0 2.0 
100.0 100.0

xsecs 1 # N Defined Materials   
#  SIGMAT---- SIGMAS---- NUSIGF----
0  1.0        0.9        0.25 

// @todo The boundary conditions should be moved to this
// file, because they belong to the problem description,
// not part of the solver strategy.

// @todo The lenght should be identified with the xsecs_id
// or node position?? 



