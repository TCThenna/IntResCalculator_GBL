########################################################
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################
# project name
PROJECT( GBLIntResCalculator )

# project version
SET( GBLINTRESCALCULATOR_VERSION_MAJOR 0 )
SET( GBLINTRESCALCULATOR_VERSION_MINOR 1 )
SET( GBLINTRESCALCULATOR_VERSION_PATCH 0 )

# some macros are redefined to keep compatability with CMake 2.6
include(${PROJECT_SOURCE_DIR}/cmake/CMakeCompatibility.cmake)

# Determine platform- and compiler-specific settings
include(${PROJECT_SOURCE_DIR}/cmake/Platform.cmake)

# Check if this is a source tarball build
IF(NOT IS_DIRECTORY ${CMAKE_SOURCE_DIR}/.git)
  SET(SOURCE_PACKAGE 1)
ENDIF(NOT IS_DIRECTORY ${CMAKE_SOURCE_DIR}/.git)

# Set package version
IF(NOT SOURCE_PACKAGE)
  # Get the version from last git tag plus numer of additional commits:
  FIND_PACKAGE(Git QUIET)
  IF(GIT_FOUND)
    EXECUTE_PROCESS(COMMAND ${GIT_EXECUTABLE} describe --tags HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} OUTPUT_VARIABLE GBLINTRESCALCULATOR_LIB_VERSION ERROR_VARIABLE GIT_ERROR)
    IF(GIT_ERROR) # git execution resulted in error message -- cannot trust result but construct version from hard-set variable
      MESSAGE(WARNING "Could not determine GBL-IntResCalc version from git -- git execution resulted in error: ${GIT_ERROR}")
      SET(GBLINTRESCALCULATOR_LIB_VERSION ${GBLINTRESCALCULATOR_VERSION_MAJOR}.${GBLINTRESCALCULATOR_VERSION_MINOR}.${GBLINTRESCALCULATOR_VERSION_PATCH})
    ELSE(GIT_ERROR)
      STRING(REGEX REPLACE "([v0-9.]+[-a-zA-Z0-9]+)-([0-9]+)-([A-Za-z0-9]+)" "\\1+\\2~\\3" GBLINTRESCALCULATOR_LIB_VERSION ${GBLINTRESCALCULATOR_LIB_VERSION})
      STRING(REGEX REPLACE "\n" "" GBLINTRESCALCULATOR_LIB_VERSION ${GBLINTRESCALCULATOR_LIB_VERSION}) # remove potential line breaks
      STRING(REGEX REPLACE "v([v0-9.]+)(.*)" "\\1" GBLINTRESCALCULATOR_LIB_BASEVERSION ${GBLINTRESCALCULATOR_LIB_VERSION})
      IF("${GBLINTRESCALCULATOR_LIB_BASEVERSION}" VERSION_LESS "${GBLINTRESCALCULATOR_VERSION_MAJOR}.${GBLINTRESCALCULATOR_VERSION_MINOR}.${GBLINTRESCALCULATOR_VERSION_PATCH}")
	MESSAGE(WARNING "You do not seem to have fetched the latest tags in your git repository -- please consider running 'git fetch' for the main remote repository.")
      ENDIF("${GBLINTRESCALCULATOR_LIB_BASEVERSION}" VERSION_LESS "${GBLINTRESCALCULATOR_VERSION_MAJOR}.${GBLINTRESCALCULATOR_VERSION_MINOR}.${GBLINTRESCALCULATOR_VERSION_PATCH}")
      EXEC_PROGRAM(git ARGS status --porcelain ${CMAKE_CURRENT_SOURCE_DIR}/src OUTPUT_VARIABLE GBLINTRESCALCULATOR_MAIN_STATUS)
      IF(GBLINTRESCALCULATOR_MAIN_STATUS STREQUAL "")
	MESSAGE(STATUS "Main GBLIntResCalculator library directory is clean (all changes have been commited).")
      ELSE(GBLINTRESCALCULATOR_MAIN_STATUS STREQUAL "")
	MESSAGE(STATUS "Main GBLIntResCalculator library directory is dirty (uncommitted changes present):\n ${GBLINTRESCALCULATOR_MAIN_STATUS}.")
	SET(GBLINTRESCALCULATOR_LIB_VERSION ${GBLINTRESCALCULATOR_LIB_VERSION}*)
      ENDIF(GBLINTRESCALCULATOR_MAIN_STATUS STREQUAL "")
    ENDIF(GIT_ERROR)
  ELSE(GIT_FOUND)
    # If we don't have git we take the hard-set version.
    SET(GBLINTRESCALCULATOR_LIB_VERSION ${GBLINTRESCALCULATOR_VERSION_MAJOR}.${GBLINTRESCALCULATOR_VERSION_MINOR}.${GBLINTRESCALCULATOR_VERSION_PATCH})
  ENDIF(GIT_FOUND)
ELSE(NOT SOURCE_PACKAGE)
  # If we don't have git we take the hard-set version.
  SET(GBLINTRESCALCULATOR_LIB_VERSION "v${GBLINTRESCALCULATOR_VERSION_MAJOR}.${GBLINTRESCALCULATOR_VERSION_MINOR}.${GBLINTRESCALCULATOR_VERSION_PATCH}")
ENDIF(NOT SOURCE_PACKAGE)
CONFIGURE_FILE("${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.cmake.h" "${CMAKE_CURRENT_BINARY_DIR}/config.h" @ONLY)
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_BINARY_DIR}")
MESSAGE("-- Determined GBLIntResCalculator version ${GBLINTRESCALCULATOR_LIB_VERSION}")

# set the path to which we will install later: default project home, can be changed using
# cmake -DINSTALL_PREFIX=/some/dir ..
SET(INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "Prefix prepended to install directories")
SET(CMAKE_INSTALL_PREFIX "${INSTALL_PREFIX}" CACHE INTERNAL "Prefix prepended to install directories" FORCE)

# additional packages to be searched for by cmake
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake )

# SET THE RPATH SO THAT OUR EXECUTABLES FIND THE TCTANALYSIS LIBRARY EVEN WHEN INSTALLED INTO NON-DEFAULT LOCATIONS
# see http://www.cmake.org/Wiki/CMake_RPATH_handling
# enable @rpath in the install name for any shared library being built
# note: it is planned that a future version of CMake will enable this by default
set(CMAKE_MACOSX_RPATH 1)
# use, i.e. don't skip the full RPATH for the build tree
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
# when building, don't use the install RPATH already
# (but later on when installing)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
# *ALTERNATIVELY*: USE RELATIVE PATHS
# see http://www.semipol.de/archives/356
# for all binaries created in a CMake project:
#SET(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/")
# for certain targets
#SET_TARGET_PROPERTIES(target 1 target2 ...
#    PROPERTIES INSTALL_RPATH "$ORIGIN/../lib:$ORIGIN/")
# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# the RPATH to be used when installing, but only if it's not a system directory
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
   SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# If the LeCroyReader package is required 
#FIND_PACKAGE( LCR REQUIRED )
#INCLUDE_DIRECTORIES( ${LCR_INCLUDE_DIR} )
#MESSAGE("--- LCR include dirs are ${LCR_INCLUDE_DIR}")
#MESSAGE("--- LCR libs are ${LCR_LIBRARIES}")
#ADD_DEFINITIONS(-DLCR_FOUND)

# Should we ever need ROOT uncomment here...
FIND_PACKAGE( ROOT REQUIRED )
INCLUDE_DIRECTORIES( ${ROOT_INCLUDE_DIR} )
MESSAGE("--- root include dirs are ${ROOT_INCLUDE_DIR}")
MESSAGE("--- root libs are ${ROOT_LIBRARIES}")
#OPTIONAL: add preprocessor defines to be able to check for ROOT in the code
ADD_DEFINITIONS(-DROOT_FOUND)

# Check for GBL
FIND_PACKAGE ( GBL REQUIRED )
INCLUDE_DIRECTORIES( ${GBL_INCLUDE_DIR} )
MESSAGE("--- GBL include dirs are ${GBL_INCLUDE_DIR}")
MESSAGE("--- GBL libs are ${GBL_LIBRARY}")
ADD_DEFINITIONS(-DGBL_FOUND)

# include directories
INCLUDE_DIRECTORIES( ./include)

# CORE LIBRARY:
SET(LIB_SOURCE_FILES
  "src/AnaTel.cc"
  )

# LCR LIBRARY:
#SET(LCR_LIB_SOURCE_FILES
#  "external/LeCroyConverter/src/LeCroy.cpp"
#  )
#ADD_LIBRARY(${PROJECT_NAME} SHARED ${LIB_SOURCE_FILES} ${LCR_LIB_SOURCE_FILES})
ADD_LIBRARY(${PROJECT_NAME} SHARED ${LIB_SOURCE_FILES} )

MESSAGE("This is the GBL_LIBRARY: ${GBL_LIBRARIES}")
# build main GBLIntResCalculator executable
add_executable(gbl-calc  src/main.cc )
target_link_libraries(gbl-calc ${PROJECT_NAME} ${ROOT_LIBRARIES} ${GBL_LIBRARIES})
# now set up the installation directory
INSTALL(TARGETS ${PROJECT_NAME} gbl-calc
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)



#############################################
# Doxygen target to generate HTML reference #
#############################################
#FIND_PACKAGE(Doxygen)
#IF(DOXYGEN_FOUND)
#  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/Doxyfile @ONLY)
#  ADD_CUSTOM_TARGET(doxygen
#    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doc/doxygen/Doxyfile
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen
#    COMMENT "Generating HTML documentation with Doxygen" VERBATIM
#    )
#ENDIF(DOXYGEN_FOUND)



#######################################################
# Print out system information useful for bug reports #
############################################### #######
include(PrintSystemInformation)


