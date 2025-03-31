set(MODULE_NAME UBWireTool)

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect Headers
set( HEADERS
  WireData.h  
  UBWireTool.h
)

# Define library and add sources
add_library( ${LIBNAME} SHARED
  WireData.cxx    
  UBWireTool.cxx
)

set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}  
  PRIVATE
)

target_link_libraries( ${LIBNAME} PUBLIC
  ${LARLITE_LIBS}
  LArCVCoreBase
  LArCVCoreDataFormat
  )

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
