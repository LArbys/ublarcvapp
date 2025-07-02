set(MODULE_NAME ContourTools)

set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect the headers
set( HEADERS
  ContourShapeMeta.h
  ContourCluster.h
  ContourClusterAlgo.h  
  #ContourAStarClusterAlgo.h
)

# Collect the sources
add_library( ${LIBNAME} SHARED
  ContourShapeMeta.cxx
  ContourCluster.cxx
  ContourClusterAlgo.cxx  
  #ContourAStarClusterAlgo.cxx
)

set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}
  ${OpenCV_INCLUDE_DIRS}
  ${GEO2D_INC_DIR}
  ${LAROPENCV_INC_DIR}
)
#target_include_directories(${LIBNAME} PUBLIC ${LARCV_INCLUDE_DIR} ${LARLITE_INC_DIRS})
target_link_libraries( ${LIBNAME} PUBLIC
  opencv_imgcodecs
  opencv_core  
  LArCVCoreBase
  LArCVCoreDataFormat
  LArCVCoreCVUtil
  ${GEO2D_LIBS}  
  ${LAROPENCV_LIBS}  
  #LArCVApp_LArOpenCVHandle
  #LArCVApp_Reco3D  
)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}
  ${OpenCV_INCLUDE_DIRS}
  ${GEO2D_INC_DIR}
)

# Generate the dictionary
ROOT_GENERATE_DICTIONARY(G__${LIBNAME} ${HEADERS} LINKDEF LinkDef.h)

# add dictionary file to sources                                                                                                                                               
target_sources( ${LIBNAME} PRIVATE G__${LIBNAME}.cxx )

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
