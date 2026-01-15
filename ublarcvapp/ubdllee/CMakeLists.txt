set(MODULE_NAME ubdllee)

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect the headers
set( HEADERS
  FixedCROIFromFlashAlgo.h
  FixedCROIFromFlashConfig.h
  dwall.h
  NuAnaMC.h
  DLInteraction.h
  MergeDLInteraction.h
  DLMerger.h
)

# Add sources
add_library( ${LIBNAME} SHARED
  dwall.cxx
  FixedCROIFromFlashAlgo.cxx
  FixedCROIFromFlashConfig.cxx
  NuAnaMC.cxx
  DLInteraction.cxx
  MergeDLInteraction.cxx
  DLMerger.cxx
)

# larlite libraries
set(LARLITE_LIBS_USED ${LARLITE_LIBS} )
set(UBLARCVAPP_LIBS_USED LArCVApp_LArliteHandler)

# Generate the shared library from the sources
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  ${LARCV_INCLUDE_DIR}
  ${LARLITE_INC_DIRS}
  ${LARCV_JSON_INCLUDE_DIR}
  PRIVATE
  ${PROJECT_SOURCE_DIR}
)  

target_link_libraries(${LIBNAME} PUBLIC
  ${LARLITE_LIBS_USED}   
  LArCVCoreBase
  LArCVCoreDataFormat
  LArCVCoreProcessor
  ${UBLARCVAPP_LIBS_USED}
)

# includes for ROOT dictionary gen
include_directories(  ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${LARCV_INCLUDE_DIR}
  ${LARCV_JSON_INCLUDE_DIR}
  ${LARLITE_INC_DIRS}
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
