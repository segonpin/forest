##
# mylib/CMakeLists.txt
##

# -------------------------------------------------------------------------
# SET THE MAIN FOLDER FOR THE SOURCES

SET(MY_DIR "./")
SET(SOURCE_DIR ${MY_DIR}/src)

# -------------------------------------------------------------------------
# SOURCES FOR THE ANGLE

SET(ANGLE_DIR ${SOURCE_DIR}/angle)
SET(ANGLE_SRC
  ${ANGLE_DIR}/quadrature_base.cc
  ${ANGLE_DIR}/quadrature_factory.cc
  ${ANGLE_DIR}/quadrature_levelsymmetric_type1.cc
  ${ANGLE_DIR}/quadrature_levelsymmetric_type2.cc
  ${ANGLE_DIR}/quadrature_gauss_chebyshev.cc
  ${ANGLE_DIR}/quadrature_chebyshev_legendre.cc
  ${ANGLE_DIR}/quadrature_gausslegendre.cc
  ${ANGLE_DIR}/moment_indexer.cc
  ${ANGLE_DIR}/spherical_harmonics.cc
  ${ANGLE_DIR}/moments_to_directions.cc
  ${ANGLE_DIR}/directions_to_moments.cc
  )

# -------------------------------------------------------------------------
# SOURCES FOR THE INPUT

SET(INPUT_DIR ${SOURCE_DIR}/input)
SET(INPUT_SRC
  ${INPUT_DIR}/input_geom.cc
  ${INPUT_DIR}/input_mat.cc
  ${INPUT_DIR}/input_settings.cc
  ${INPUT_DIR}/input.cc
  )

# -------------------------------------------------------------------------
# SOURCES FOR THE GEOMETRY

SET(GEOMETRY_DIR ${SOURCE_DIR}/geometry)
SET(GEOMETRY_SRC
  ${GEOMETRY_DIR}/real_pin.cc
  ${GEOMETRY_DIR}/real_assembly.cc
  ${GEOMETRY_DIR}/real_core.cc
  ${GEOMETRY_DIR}/prob_geom.cc
  ${GEOMETRY_DIR}/path_to_gridgenerator_dealii.cc
  )

# -------------------------------------------------------------------------
# SOURCES FOR THE NEUTRONICS

SET(NEUTRONICS_DIR ${SOURCE_DIR}/neutronics)
SET(NEUTRONICS_SRC
  ${NEUTRONICS_DIR}/scattering_source.cc
  ${NEUTRONICS_DIR}/scattering_factory.cc
  ${NEUTRONICS_DIR}/scattering_matrix_built.cc
  ${NEUTRONICS_DIR}/scattering_matrix_free.cc
  ${NEUTRONICS_DIR}/fission_source_matrix_built.cc
  ${NEUTRONICS_DIR}/fission_source_matrix_free.cc
  ${NEUTRONICS_DIR}/fission_source_factory.cc
  ${NEUTRONICS_DIR}/fission_spectrum_matrix_built.cc
  ${NEUTRONICS_DIR}/fission_spectrum_matrix_free.cc
  ${NEUTRONICS_DIR}/fission_spectrum_factory.cc
  ${NEUTRONICS_DIR}/diff_factory.cc
  ${NEUTRONICS_DIR}/diff_op.cc
  ${NEUTRONICS_DIR}/diff_op_matrix_built.cc
  ${NEUTRONICS_DIR}/diff_op_matrix_free.cc
  ${NEUTRONICS_DIR}/trans_factory.cc
  ${NEUTRONICS_DIR}/trans_op.cc
  ${NEUTRONICS_DIR}/trans_op_matrix_built.cc
  ${NEUTRONICS_DIR}/trans_op_matrix_free.cc
  ${NEUTRONICS_DIR}/trans_op_mf_pre.cc
  ${NEUTRONICS_DIR}/trans_bc.cc
  ${NEUTRONICS_DIR}/boundary_values_factory.cc
  ${NEUTRONICS_DIR}/boundaryvalues.cc
  ${NEUTRONICS_DIR}/downstream_new.cc
  ${NEUTRONICS_DIR}/neutronicmodule.cc
  ${NEUTRONICS_DIR}/state.cc
  ${NEUTRONICS_DIR}/extract_bv.cc
  ${NEUTRONICS_DIR}/path_to_intergridmap_dealii.cc
  ${NEUTRONICS_DIR}/homogenization.cc
  ${NEUTRONICS_DIR}/hom_data.cc
  #${NEUTRONICS_DIR}/nonlinearresidual.cc
  #${NEUTRONICS_DIR}/jfnk.cc
  )

# -------------------------------------------------------------------------
# SOURCES FOR THE UTILS

SET(UTILS_DIR ${SOURCE_DIR}/utils)
SET(UTILS_SRC
  ${UTILS_DIR}/forest_utils_base.cc
  ${UTILS_DIR}/forest_utils_dealii.cc
  ${UTILS_DIR}/forest_utils_timer.cc
  ${UTILS_DIR}/forest_utils_memory.cc
  ${UTILS_DIR}/cmdoption.cc
  )

# -------------------------------------------------------------------------
# SOURCES FOR THE OUTPUT
# @TODO Change "output" for "postprocessing" and add more functionalities? 

SET(OUTPUT_DIR ${SOURCE_DIR}/output)
SET(OUTPUT_SRC
  ${OUTPUT_DIR}/forest_post_processing.cc
  )

# -------------------------------------------------------------------------
# SOURCES FOR THE NOT-CLASsIFIED-YET
# @TODO classify this files in the right folder. 

SET(ALGEBRA_DIR ${SOURCE_DIR}/algebra)
SET(ALGEBRA_SRC
  ${ALGEBRA_DIR}/forest_vector.cc
  ${ALGEBRA_DIR}/sn_vector.cc
  ${ALGEBRA_DIR}/eq_factory.cc
  ${ALGEBRA_DIR}/eq_transport.cc
  ${ALGEBRA_DIR}/eq_diffusion.cc
  ${ALGEBRA_DIR}/solver_eq.cc
  ${ALGEBRA_DIR}/eigen_solver.cc
  ${ALGEBRA_DIR}/eigen_prob_factory.cc
  ${ALGEBRA_DIR}/eigen_standard.cc
  ${ALGEBRA_DIR}/eigen_fission_density.cc
  ${ALGEBRA_DIR}/eigen_solver_slepc.cc
  )

# -------------------------------------------------------------------------
# Declare all source files the target consists of:
# -------------------------------------------------------------------------

SET(TARGET_SRC
  ${INPUT_SRC}
  ${GEOMETRY_SRC}
  ${ANGLE_SRC}
  ${NEUTRONICS_SRC}
  ${OUTPUT_SRC}
  ${UTILS_SRC}
  ${ALGEBRA_SRC}
  )

# -------------------------------------------------------------------------

INCLUDE_DIRECTORIES(include)
ADD_LIBRARY(mylib ${TARGET_SRC})
set_target_properties(mylib PROPERTIES PREFIX "lib")
DEAL_II_SETUP_TARGET(mylib)

# -------------------------------------------------------------------------
