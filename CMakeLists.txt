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
  find_package(OpenCV REQUIRED)
  if (NOT OpenCV_FOUND)
    message("OpenCV NOT FOUND")
    set(USE_OPENCV 0)
  endif()
  message("OpenCV FOUND")
else()
  message("OpenCV not activated")
endif()

#---larbys GEO2D
# use it if env variables defined
if (DEFINED ENV{GEO2D_BASEDIR})
  set(LARCV_USE_GEO2D 1)
  
  # unless explicitly told no
  if (DEFINED ENV{USE_GEO2D} AND NOT $USE_GEO2D)
    set(LARCV_USE_GEO2D 0)
    message("DISABLING modules that need Geo2D")
  else()
    file( GLOB GEO2D_LIBS "$ENV{GEO2D_LIBDIR}/*.so" )
    set( GEO2D_INC_DIR "$ENV{GEO2D_INCDIR}" )
    set( GEO2D_LIB_DIR "$ENV{GEO2D_LIBDIR}" )
  endif()
else()
  set(LARCV_USE_GEO2D 0)
  message("DISABLING modules needing Geo2D")  
endif()


#---LArLite
set( LARLITE_INC_DIRS "$ENV{LARLITE_COREDIR}" "$ENV{LARLITE_BASEDIR}" "$ENV{LARLITE_USERDEVDIR}" )
set( LARLITE_LIB_DIR  "$ENV{LARLITE_LIBDIR}" )
set( HAS_LARLITE 1 )
file( GLOB LARLITE_LIBS "$ENV{LARLITE_LIBDIR}/libLArLite*.so" )

#---LArOpenCV
set( LAROPENCV_INC_DIR "$ENV{LAROPENCV_BASEDIR}" )
set( LAROPENCV_LIB_DIR "$ENV{LARLITE_LIBDIR}" )
file( GLOB LAROPENCV_LIBS "$ENV{LARLITE_LIBDIR}/libLArOpenCV_*.so" )

#---Cilantro
set( cilantro_DIR $ENV{CILANTRO_LIB_DIR} )
find_package( cilantro REQUIRED )


## Setup

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

#---Define Variables for subdirectory cmakelists
set(UBLARCVAPP_BASE_DIR "${CMAKE_CURRENT_SOURCE_DIR}")

## Add Source

#---Add base of source
add_subdirectory(ublarcvapp)

#---Configuration export
include(CMakePackageConfigHelpers)
set(UBLARCVAPP_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}/include")
set(UBLARCVAPP_LIBRARY_DIR "${CMAKE_INSTALL_PREFIX}/lib")
set(TARGETS_EXPORT_NAME "ublarcvappTargets")
configure_package_config_file("cmake/Config.cmake.in"
  "${CMAKE_BINARY_DIR}/cmake/ublarcvappConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ublarcvapp"
  PATH_VARS UBLARCVAPP_INCLUDE_DIR UBLARCVAPP_LIBRARY_DIR UBLARCVAPP_BASE_DIR)
write_basic_package_version_file("${CMAKE_BINARY_DIR}/cmake/ublarcvappConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion)

## Export Targets

#---Create info to include LArCV in other cmake files
install(EXPORT ublarcvapp
        FILE ublarcvappTargets.cmake
        DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ublarcvapp")

#---Create info to include LArCV in other cmake files
install(FILES "${CMAKE_BINARY_DIR}/cmake/ublarcvappConfig.cmake"
  "${CMAKE_BINARY_DIR}/cmake/ublarcvappConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/ublarcvapp")







