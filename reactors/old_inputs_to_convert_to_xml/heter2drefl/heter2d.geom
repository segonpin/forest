# ##########################################################
# Heterogeneous square
# The reference results have been obtained with dragon, usging
# degree three (d:=3), no refinement (s:=1) and type 2
# quadrature set (QUAD = 2);
# The reference for sn = 02 is keff = 5.755775E-01
# The reference for sn = 04 is keff = 5.806699E-01
# The reference for sn = 06 is keff = 5.900330E-01
# The reference for sn = 08 is keff = 5.941040E-01
# The reference for sn = 12 is keff = 5.957445E-01
# The reference for sn = 16 is keff = 5.963024E-01
# The reference for sn = 20 is keff = 5.964814E-01
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
n_nodes 4 4
n_xsecs 2
n_groups 1

# material distribution
Materials
1   1   1   1
1   0   1   0
1   1   1   1
1   0   1   0

# I have to change that to have
# step_sizes
# x_lenght 2.7 2.4 2.7 2.4 2.7
# y_lenght 2.7 2.4 2.7 2.4 2.7
# z_lenght 1.0 1.0 1.0 1.0 1.0

step_length 
2.7 2.4 2.7 1.2
2.7 2.4 2.7 1.2

# #
# The cross sections description follows below for different mixtures
xsecs 2 # N Defined Materials
# #
# the first mixture
MIX 0
SIGMAT 0.41666700
SIGMAS 0.33400000
NUSIGF 0.17800000
CHI    1.00000000
# #
# now the second mixture
MIX 1
SIGMAT 0.37037000
SIGMAS 0.33400000
NUSIGF 0.00000000
CHI    1.00000000

# #
# 0 means vacuum, 1 means incoming flux that should be specified
# by face (more complicated definition of the incoming flux should be
# available, like projection of something else), and 2 means 
# specular reflection boundary condition.
# the order is Left, Right, Front, Back, Bottom and Top face, (-x, x, -y, y, -z, z)
Boundary_Conditions
bc_id 0 2 0 2
in_id 0 0 0 0





