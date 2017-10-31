# ##########################################################
# Heterogeneous slab
# The reference results have been obtained with dragon,
# using degree three (d:=3),  and no refinement (s:=1)
# second type quadrature (QUAD 2):
# The reference for sn = 02 is keff = 1.105293
# The reference for sn = 04 is keff = 1.147933
# The reference for sn = 10 is keff = 1.161048
# The reference for sn = 32 is keff = 1.162142
# The reference for sn = 64 is keff = 1.162215
# The reference for sn = 96 is keff = 1.16224
# This problem is provided in the following paper:
#
# M. Capilla,  CF. Talavera,  D. Ginestar,  G. Verdu, 
# "A nodal collocation method for the calculation of
# the lambda modes of the P_L equations";
# Annals of Nuclear Energy, 32,  (2005), 1825--1853.
# ##########################################################


# Definitions of materials

n_nodes 7 
n_xsecs 2 

Materials
 1   0   1   0   1   0   1

step_length
2.7 2.4 2.7 2.4 2.7 2.4 2.7

xsecs 2 # N Defined Materials   
#  SIGMAT---- SIGMAS---- NUSIGF----
MIX 0
SIGMAT 0.41666700
SIGMAS 0.33400000
NUSIGF 0.17800000
CHI    1.00000000
# now the second mixture
MIX 1
SIGMAT 0.37037000
SIGMAS 0.33400000
NUSIGF 0.00000000
CHI    1.00000000

// @todo The boundary conditions should be moved to this
// file, because they belong to the problem description,
// not part of the solver strategy.

// @todo The lenght should be identified with the xsecs_id
// or node position?? 



