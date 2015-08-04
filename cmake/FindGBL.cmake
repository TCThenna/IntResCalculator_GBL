# - Find GBL installation
# This module tries to find the libGBL installation on your system.
# Once done this will define
#
#  GBL_FOUND - system has GBL installed
#  GBL_INCLUDE_DIR - ~ the GBL include directory 
#  GBL_LIBRARY - Link these to use GBL

MESSAGE(STATUS "Looking for GBL...")

FIND_PATH(GBL_INCLUDE_DIR GblPoint.h
  HINTS "${PROJECT_SOURCE_DIR}/include"
)

FIND_LIBRARY(GBL_LIBRARY NAMES libGBL.so
  HINTS "${PROJECT_SOURCE_DIR}/lib"
)

set(GBL_LIBRARIES ${GBL_LIBRARY})

IF (GBL_LIBRARY)
  IF(GBL_INCLUDE_DIR)
    set(GBL_FOUND TRUE)
    set(GBL_LIBRARIES "${GBL_LIBRARIES}/libGbl.so")
    MESSAGE(STATUS "Found libGBL: ${GBL_INCLUDE_DIR}, ${GBL_LIBRARIES}")
  ELSE(GBL_INCLUDE_DIR)
    set(GBL_FOUND FALSE)
    MESSAGE(STATUS "GBL headers NOT FOUND. Make sure to install the headers.")
  ENDIF(GBL_INCLUDE_DIR)
ELSE (GBL_LIBRARY)
    set(GBL_FOUND FALSE)
    MESSAGE(STATUS "libGBL NOT FOUND.")
ENDIF (GBL_LIBRARY)

set(GBL_INCLUDE_DIR
    ${GBL_INCLUDE_DIR}
)
