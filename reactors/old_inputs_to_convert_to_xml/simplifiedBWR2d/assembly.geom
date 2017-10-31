# ##########################################################
# Heterogeneous square
# The reference results have been obtained with dragon, usging
# degree three (d:=?), no refinement (s:=?) and type ? 
# quadrature set (QUAD = ? );
# The reference for sn = 02 is keff = ? 
# The reference for sn = 04 is keff = ? 
# The reference for sn = 06 is keff = ? 
# The reference for sn = 08 is keff = ? 
# The reference for sn = 12 is keff = ? 
# The reference for sn = 16 is keff = ? 
# The reference for sn = 20 is keff = ? 
#
# This problem is provided by:
#
# Roberts, Jeremy A., and Benoit Forget.
# "Multigroup diffusion preconditioners for multiplying 
# fixed-source transport problems."
# Journal of Computational Physics (2014).
# ##########################################################


# Definitions of materials

n_nodes 19 19
n_xsecs 2 
n_groups 3

Materials
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 0 1 1 1 1 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 1 1 1 1 1 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 1 1 1 1 0 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1 0 1
1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1


# I have to change that to have
# step_sizes
# x_lenght 2.7 2.4 2.7 2.4 2.7
# y_lenght 2.7 2.4 2.7 2.4 2.7
# z_lenght 1.0 1.0 1.0 1.0 1.0

step_length 
0.825 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.825 
0.825 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.25 1.150 0.825 

xsecs 2 # N Defined Materials   
#  SIGMAT---- SIGMAS---- NUSIGF----
# Mixture number zero (the fuel)
MIX 0
SIGMAT 0.2822058997  0.4997685502  0.4323754911
# SIGMAF 0.0028231045  0.0096261203  0.1123513981
# NU     2.7202775245  2.4338148428  2.4338000000
CHI    0.9996892490  0.0003391680  0.0000000000
SIGMAS 0.2760152893  0.0000000000  0.0000000000
       0.0011230014  0.4533430274  0.0000378305
       0.0000000000  0.0014582502  0.2823864370
NUSIGF 0.0076796277  0.0234281944  0.2734408326
# Mixture number one (the moderator)
MIX 1
SIGMAT 0.3333694980  0.5887110656  1.6512246291
SIGMAS 0.2432946408  0.0000000000  0.0000000000
       0.0898364840  0.4413753398  0.0000203109
       0.0000387911  0.1465683257  1.6300848232
NUSIGF 0.0000000000  0.0000000000  0.0000000000
CHI    0.0000000000  0.0000000000  0.0000000000



# #
# 0 means vacuum, 1 means incoming flux that should be specified
# by face (more complicated definition of the incoming flux should be
# available, like projection of something else), and 2 means 
# specular reflection boundary condition.
# the order is Left, Right, Front, Back, Bottom and Top face, (-x, x, -y, y, -z, z)
Boundary_Conditions
bc_id 0 0 0 0
in_id 0 0 0 0





