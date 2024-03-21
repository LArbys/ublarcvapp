set(MODULE_NAME MCTools)

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect the headers
set( HEADERS MCPos2ImageUtils.h
  MCPixelPGraph.h
  TruthTrackSCE.h
  TruthShowerTrunkSCE.h
  crossingPointsAnaMethods.h
  NeutrinoVertex.h
  LArbysMC.h
  NeutrinoPixelFilter.h
  FlashMatcher.h
  FlashMatcherV2.h  
  MCPixelPMap.h )

# Add sources
add_library( ${LIBNAME} SHARED MCPos2ImageUtils.cxx
  MCPixelPGraph.cxx
  TruthTrackSCE.cxx
  TruthShowerTrunkSCE.cxx
  crossingPointsAnaMethods.cxx
  NeutrinoVertex.cxx
  LArbysMC.cxx
  NeutrinoPixelFilter.cxx
  FlashMatcher.cxx
  FlashMatcherV2.cxx  
  MCPixelPMap.cxx )

# larlite libraries
set(LARLITE_LIBS_USED ${LARLITE_LIBS} )

# includes for ROOT dictionary gen
include_directories(${LARCV_INCLUDE_DIR} ${LARCV_JSON_INCLUDE_DIR} ${LARLITE_INC_DIRS} ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/ublarcvapp/${MODULE_NAME})

if (HAS_LARLITE)
  set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -DHAS_LARLITE")
  #message("build ubimagemod with larlite")
else()
  #message("build ubimagemod without larlite")
endif()

# Generate the dictionary
ROOT_GENERATE_DICTIONARY(G__${LIBNAME} ${HEADERS} LINKDEF LinkDef.h)
target_sources( ${LIBNAME} PRIVATE G__${LIBNAME}.cxx )

# Generate the shared library from the sources
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_include_directories(${LIBNAME} PUBLIC ${LARCV_INCLUDE_DIR} ${LARLITE_INC_DIRS} ${LARCV_JSON_INCLUDE_DIR} )
target_link_libraries(${LIBNAME} LArCVCoreBase LArCVCoreDataFormat LArCVCoreProcessor ${LARLITE_LIBS_USED})

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
