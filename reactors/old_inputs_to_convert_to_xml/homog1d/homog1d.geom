# ##########################################################
# Homogeneous slab 10.0 cm (it can be a 2.0 cm or 10.0 cm slab)
# The reference for sn = 96 is keff = 0.662951 (2cm slab)
# The reference for sn = 96 is keff = 2.00193 (10cm slab)
# This problem is provided in the following paper:
#
# M. Capilla,  CF. Talavera,  D. Ginestar,  G. Verdu, 
# "A nodal collocation method for the calculation of
# the lambda modes of the P_L equations";
# Annals of Nuclear Energy, 32,  (2005), 1825--1853.
# ##########################################################

# Definitions of materials
n_nodes 5
n_xsecs 1

Materials
 0   0   0   0   0 

step_length
2.0 2.0 2.0 2.0 2.0 
# 0.4 0.4 0.4 0.4 0.4 

xsecs 1 # Defined Materials   
#  SIGMAT---- SIGMAS---- NUSIGF----
0  1.0        0.9        0.25


// @todo The boundary conditions should be moved to this
// file, because they belong to the problem description,
// not part of the solver strategy.

// @todo The lenght should be identified with the xsecs_id
// or node position?? 


