cmake_minimum_required(VERSION 3.10)
project(ublarcvapp VERSION 1.1.0)

#---Install location by default is in the build directory
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}" CACHE PATH "${CMAKE_BUILD_DIR}" FORCE )
endif()

# CFLAGS
SET(CMAKE_C_FLAGS "-std=gnu99 -pedantic -fPIC -pedantic -fopenmp -mtune=native")
SET(CMAKE_CXX_FLAGS "-Wall -fPIC -pedantic -fopenmp -mtune=native")

## Dependencies

#---CERN ROOT (required)
find_package(ROOT CONFIG REQUIRED)
# we set the compiler version based on the one used for ROOT
execute_process(
    COMMAND bash "-c" "root-config --features | awk '{print $1}'"
    OUTPUT_VARIABLE ROOT_CXX_STANDARD )
if (${ROOT_CXX_STANDARD} MATCHES "cxx17")
  message( STATUS "using cxx17")
  set(CMAKE_CXX_STANDARD 17)
elseif (${ROOT_CXX_STANDARD} MATCHES "cxx14")
  message( STATUS "using cxx14")
  set(CMAKE_CXX_STANDARD 14)
else ()
  message( STATUS "using cxx11")  
  set(CMAKE_CXX_STANDARD 11)
endif()
    
#---LArCV (required)
set( larcv_DIR "$ENV{LARCV_LIBDIR}/cmake/larcv")
find_package( larcv REQUIRED )
#   -- json header included in larcv
set(LARCV_JSON_INCLUDE_DIR "$ENV{LARCV_BASEDIR}/extern/json/single_include")

#---OPENCV (optional)
if (USE_OPENCV)
  set(OpenCV_DIR $ENV{OPENCV_LIBDIR}/../share/OpenCV)
  find_package(OpenCV REQUIRED)
  if (NOT OpenCV_FOUND)
    message("OpenCV NOT FOUND")
    set(USE_OPENCV 0)
  endif()
  #message("OpenCV FOUND")
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

#---Eigen3
find_package(Eigen3 REQUIRED)

#---LArLite
set( LARLITE_INC_DIRS "$ENV{LARLITE_INCDIR}" )
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

#---Boost
find_package( Boost 1.65 COMPONENTS graph program_options )
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON) 
set(Boost_USE_STATIC_RUNTIME OFF)
if (Boost_FOUND)
    include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
    message(STATUS "Using boost include dirs: ${Boost_INCLUDE_DIRS}")
    message(STATUS "Using boost library dirs: ${Boost_LIBRARY_DIRS}")
    message(STATUS "Using boost libraries: ${Boost_LIBRARIES}")
    message(STATUS "Using boost thread library: ${Boost_THREAD_LIBRARY}")
    link_directories(${Boost_LIBRARY_DIRS})
endif (Boost_FOUND)

#INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

#--- glog
#find_package(glog REQUIRED)


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







