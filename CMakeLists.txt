cmake_minimum_required(VERSION 3.10)
project(ublarcvapp VERSION 1.1.0)

#---Install location by default is in the build directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "${CMAKE_BUILD_DIR}" FORCE )
endif()

## Dependencies

#---CERN ROOT (required)
find_package(ROOT CONFIG REQUIRED)

#---LArCV (required)
set( larcv_DIR "$ENV{LARCV_LIBDIR}/cmake/larcv")
find_package( larcv REQUIRED )

#---OPENCV (optional)
if (USE_OPENCV)
  find_package(OpenCV)
  if (NOT OpenCV_FOUND)
    set(USE_OPENCV 0)
  endif()
endif()

#---larbys GEO2D
# use it if env variables defined
if (DEFINED ENV{GEO2D_BASEDIR})
  set(LARCV_USE_GEO2D 1)
  file( GLOB GEO2D_LIBS "$ENV{GEO2D_LIBDIR}/*.so" )
  set( GEO2D_INC_DIR "$ENV{GEO2D_INCDIR}" )
  set( GEO2D_LIB_DIR "$ENV{GEO2D_LIBDIR}" )
else()
  set(LARCV_USE_GEO2D 0)
endif()

# unless explicitly told no
if (NOT USE_GEO2D)
  set(LARCV_USE_GEO2D 0)
endif()

#---LArLite
set( LARLITE_INC_DIRS "$ENV{LARLITE_COREDIR}" "$ENV{LARLITE_BASEDIR}" "$ENV{LARLITE_USERDEVDIR}" )
set( LARLITE_LIB_DIR  "$ENV{LARLITE_LIBDIR}" )
file( GLOB LARLITE_LIBS "$ENV{LARLITE_LIBDIR}/libLArLite*.so" )

#---LArOpenCV
set( LAROPENCV_INC_DIR "$ENV{LAROPENCV_BASEDIR}" )
set( LAROPENCV_LIB_DIR "$ENV{LARLITE_LIBDIR}" )
file( GLOB LAROPENCV_LIBS "$ENV{LARLITE_LIBDIR}/libLArOpenCV_*.so" )


## Setup

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---Define Variables for subdirectory cmakelists
set(UBLARCVAPP_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

## Add Source

#---Add base of source
add_subdirectory(ublarcvapp)

## Export Targets

#---Create info to include LArCV in other cmake files
install(EXPORT ublarcvapp
        FILE ublarcvappTargets.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ublarcvapp")






