##
# testing/CMakeLists.txt
##
# -------------------------------------------------------------------------

ADD_EXECUTABLE(main.exe main.cc)
INCLUDE_DIRECTORIES(${mylib_includes} ${Boost_INCLUDE_DIRS} ${DEAL_II_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(main.exe mylib)
# we add this to link with lapack and blas
TARGET_LINK_LIBRARIES(main.exe lapack blas)

DEAL_II_SETUP_TARGET(main.exe)

SET(RUNCASES 2x2pins 3x3pins c5g7_inner_uo2 c5g7_mox c5g7_outer_uo2
	     c5g7_s10 c5g7_s6 c5g7 heter1d heter2d homog1d homog2d
	     lattice2d minicore2d pin2d simplifiedBWR2d_grid
	     simplifiedBWR2d_pins GeorgiaTech1D
    )
FOREACH(elem ${RUNCASES})
  ADD_CUSTOM_TARGET(run${elem} COMMAND env CMAKE_SOURCE_DIR=${CMAKE_SOURCE_DIR}
    sh ${elem}.sh
    DEPENDS main.exe
  )
ENDFOREACH()
