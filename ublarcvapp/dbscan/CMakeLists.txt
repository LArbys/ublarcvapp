set(MODULE_NAME dbscan)

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect the headers
set( HEADERS
  DBScan.h
  DBScanTypes.h
  kdtree.h
  sDBScan.h
)

# Collect the sources
add_library( ${LIBNAME} SHARED
  DBScan.cxx
  kdtree.cxx
)

set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  PRIVATE
  ${PROJECT_SOURCE_DIR}
)
#target_include_directories(${LIBNAME} PUBLIC ${LARCV_INCLUDE_DIR} ${LARLITE_INC_DIRS})
target_link_libraries( ${LIBNAME} PUBLIC
  LArCVCoreBase
  LArCVCoreDataFormat
)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${LARCV_INCLUDE_DIR}
  ${cilantro_INCLUDE_DIRS}
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
