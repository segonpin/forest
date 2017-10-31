##
# top level CMakeLists.txt
##

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

if (${CMAKE_EXTRA_GENERATOR} MATCHES "Eclipse CDT4")
   set(CMAKE_CXX_COMPILER_ARG1 "-std=c++11" CACHE STRING "C++ version for eclipse" FORCE)
   set(CMAKE_ECLIPSE_VERSION "4.6" CACHE STRING "Eclipse version" FORCE)
endif (${CMAKE_EXTRA_GENERATOR} MATCHES "Eclipse CDT4")

# -------------------------------------------------------------------------

# This is just a message to know that the cmake is working
MESSAGE("====================================================")
MESSAGE("============ Configuring FOREST ====================")
MESSAGE("====================================================")

# -------------------------------------------------------------------------

# The minimum required version that we need of cmake because of deal.ii
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

# -------------------------------------------------------------------------

ADD_DEFINITIONS(-DPROJECT_DIRECTORY=\"${CMAKE_SOURCE_DIR}\")

# We look for deal.ii
FIND_PACKAGE(deal.II QUIET
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE("DEAL_II_DIR = ${DEAL_II_DIR}")
  MESSAGE("deal.II_DIR = ${deal.II_DIR}")
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
    set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif(CCACHE_FOUND)


# We print a message saying the version found
MESSAGE(STATUS "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'")
 
ADD_CUSTOM_TARGET(debug
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Debug ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Debug"
  )

ADD_CUSTOM_TARGET(release
  COMMAND ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=Release ${CMAKE_SOURCE_DIR}
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target all
  COMMENT "Switch CMAKE_BUILD_TYPE to Release"
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0 -Wno-unused-local-typedefs")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations") #avoid boost warning for auto_ptr

# -------------------------------------------------------------------------

# Here is the name of the project. I dont know the name yet, but without
# this line the Cmake does not work properly
PROJECT(FOREST)

# Options for shared library. On or Off.
#set(BUILD_SHARED_LIBS YES CACHE BOOL        "Enable shared libraries.")
#set(LIB_TYPE STATIC)
#if(BUILD_SHARED_LIBS)
#  set(LIB_TYPE SHARED)
#endif()

# Options for shared library. On or Off.
set(BUILD_SHARED_LIBS YES CACHE BOOL        "Enable shared libraries.")
set(LIB_TYPE STATIC)
if(BUILD_SHARED_LIBS)
  set(LIB_TYPE SHARED)
endif()


# Options for unit tests. On or Off.
set(FOREST_ENABLE_TEST FALSE CACHE BOOL        "Enable the units tests.")
# Options for Python bindings. On or Off.
set(FOREST_ENABLE_PYTHON TRUE CACHE BOOL      "Enable Python bindings via SWIG.")
# Options for PETSc. On or Off.
#set(FOREST_ENABLE_PETSC FALSE CACHE BOOL       "Build callow with PETSc.")
# Options for SLEPc. On or Off.
#set(FOREST_ENABLE_SLEPC FALSE CACHE BOOL       "Build callow SLEPc.")
# Options for Boost. On or Off.
#set(FOREST_ENABLE_BOOST YES CACHE BOOL       "Enable BOOST.")


SET(FOREST_ENABLE_BOOST ON CACHE BOOL "Enable BOOST.")
if(FOREST_ENABLE_BOOST)
  MESSAGE("====================================================")
  MESSAGE("Adding BOOST")
  set(Boost_USE_STATIC_LIBS OFF)
  set(Boost_USE_MULTITHREAD OFF)
  set(Boost_ADDITIONAL_VERSIONS "1.50" "1.50.0")
  find_package(Boost 1.50 COMPONENTS serialization)
  if(Boost_FOUND)
    MESSAGE("BOOST found")
    include_directories(${Boost_INCLUDE_DIRS})
  else(Boost_FOUND)
    set(FOREST_ENABLE_BOOST OFF)
    message(STATUS "Disabling BOOST since it was not found.")
    # Reset to avoid downstream errors due to NOTFOUND's.
    set(Boost_LIBRARIES "")
    set(Boost_INCLUDE_DIRS "")
  endif(Boost_FOUND)
  MESSAGE("====================================================")
else(FOREST_ENABLE_BOOST)
  MESSAGE("BOOST not enable")
endif()

#include(ProcessorCount)
#ProcessorCount(NPROC)
#MATH(EXPR NPROC "${NPROC}-1")
#set(CTEST_BUILD_FLAGS -j${NPROC})
#set(ctest_test_args ${ctest_test_args} PARALLEL_LEVEL ${NPROC})


# -------------------------------------------------------------------------
# I have to add the library
MESSAGE("====================================================")
MESSAGE("Adding mylib")
ADD_SUBDIRECTORY(mylib)
set(mylib_includes ${CMAKE_SOURCE_DIR}/mylib/include)
include_directories(${mylib_includes})
MESSAGE("====================================================")

# -------------------------------------------------------------------------
# I have to add all the examples to this variable
#SET(BUILD_EXAMPLES FALSE)
#IF(${BUILD_EXAMPLES})
#  SET(EXAMPLES
#    examples/test1
#    examples/test2
#    )
#  FOREACH(element ${EXAMPLES})
#    MESSAGE("Adding ${element}")
#    ADD_SUBDIRECTORY(${element})
#  ENDFOREACH()
#ENDIF()
#MESSAGE("====================================================")
# -------------------------------------------------------------------------
# I have to add the working directory
SET(BUILD_TESTING TRUE)
IF(${BUILD_TESTING})
  SET(TESTING testing)
  MESSAGE("Adding ${TESTING}")
  ADD_SUBDIRECTORY(${TESTING})
ENDIF()
MESSAGE("====================================================")


#set_property(TARGET ${TARGET} PROPERTY LINK_DEPENDS ${MYLIB})

# -------------------------------------------------------------------------
