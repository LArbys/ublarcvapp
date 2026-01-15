set(MODULE_NAME UBPhotonLib)

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect the headers
set( HEADERS
  UBPhotonLib.h
  PhotonVisibilityEstimator.h
)

# Add sources
add_library( ${LIBNAME} SHARED
  UBPhotonLib.cxx
  PhotonVisibilityEstimator.cxx
)

# larlite libraries
set(LARLITE_LIBS_USED ${LARLITE_LIBS} )

set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}  
  PRIVATE
  ${PROJECT_SOURCE_DIR}
)

target_link_libraries( ${LIBNAME} PUBLIC
  ${LARLITE_LIBS_USED}  
  LArCVCoreBase
  LArCVCoreDataFormat
  LArCVCoreProcessor
)

# includes for ROOT dictionary gen
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}
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
