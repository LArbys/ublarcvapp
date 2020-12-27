set(MODULE_NAME DLTagger)

# DEFINE THE TARGET
set( LIBNAME LArCVApp_${MODULE_NAME} )

set( HEADERS DLTagger.h
  DLTaggerProcess.h
  DLCheaterTagger.h
  DLSSNetOverlayCheater.h )

add_library( ${LIBNAME} SHARED
  DLTagger.cxx
  DLTaggerProcess.cxx
  DLCheaterTagger.cxx
  DLSSNetOverlayCheater.cxx )

# create a property to hold our headers
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

# DEFINE SOME VARS THAT MIGHT BE USEFUL FOR SUBMODULES

# larlite includes
set(LARLITE_INCLUDES $ENV{LARLITE_BASEDIR} $ENV{LARLITE_USERDEVDIR} $ENV{LARLITE_COREDIR} $ENV{LARLITE_USERDEVDIR}/BasicTool )

# larlite libraries
#set(LARLITE_LIBS_USED SelectionTool_LEEPreCuts)

# larcv libs used
set(LARCV_LIBS_USED LArCVCoreBase LArCVCoreDataFormat LArCVCoreCVUtil)

# ublarcvapp libraries
set(UBLARCVAPP_LIBS_USED LArCVApp_ContourTools LArCVApp_UBImageMod LArCVApp_Reco3D)

# sub-directories
add_subdirectory( mrcnnmatch )

message("headers: ${HEADERS}")
message("project source dir: ${PROJECT_SOURCE_DIR}")
message("eigen3: ${EIGEN3_INCLUDE_DIRS}")

get_target_property(ALL_HEADERS ${LIBNAME} PUBLIC_HEADER)

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/mrcnnmatch>
  ${LARCV_INCLUDE_DIR}
  ${LARLITE_INCLUDES}
  ${cilantro_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
  ${GEO2D_INC_DIR}
  ${Boost_INCLUDE_DIR}
  PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/mrcnnmatch
  ${PROJECT_SOURCE_DIR}
  )
target_link_libraries(${LIBNAME} ${LARLITE_LIBS_USED} ${LARCV_LIBS_USED} ${UBLARCVAPP_LIBS_USED} ${LAROPENCV_LIBS} ${GEO2D_LIBS} opencv_core)

include_directories( ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/mrcnnmatch ${LARCV_INCLUDE_DIR} ${LARLITE_INCLUDES} $ENV{CILANTRO_INC_DIR} ${EIGEN3_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS} )
ROOT_GENERATE_DICTIONARY(G__${LIBNAME} ${ALL_HEADERS} LINKDEF LinkDef.h)
target_sources( ${LIBNAME} PRIVATE G__${LIBNAME}.cxx )


message("public headers: ${ALL_HEADERS}")

# install libraries and headers
install(TARGETS ${LIBNAME}
  EXPORT ublarcvapp
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/ublarcvapp/${MODULE_NAME})

# install files needed for dictionary use in CINT/CLING/PYTHON
install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}.rootmap
  DESTINATION lib  )
