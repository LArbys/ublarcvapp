set(MODULE_NAME UBImageMod)

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect Headers
set( HEADERS
  EmptyChannelAlgo.h
  InfillDataCropper.h
  InfillImageStitcher.h
  InfillSparsifyImage.h
  TrackImageMask.h
  UBCropLArFlow.h
  UBCropMask.h
  UBSplitDetector.h
  )

# Define library and add sources
add_library( ${LIBNAME} SHARED
  EmptyChannelAlgo.cxx
  InfillDataCropper.cxx
  InfillImageStitcher.cxx
  InfillSparsifyImage.cxx
  TrackImageMask.cxx
  UBCropLArFlow.cxx
  UBCropMask.cxx
  UBSplitDetector.cxx  
  )

set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  ${LARLITE_INCLUDE_DIR}
  ${LARCV_INCLUDE_DIR}  
  PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}  
)

target_link_libraries( ${LIBNAME} PUBLIC
  ${LARLITE_LIBS}
  LArCVCoreBase
  LArCVCoreDataFormat
  LArCVCoreProcessor
)

# includes for ROOT dict generator
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${LARCV_INCLUDE_DIR}
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
