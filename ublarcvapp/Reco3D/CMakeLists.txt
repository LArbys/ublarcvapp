set(MODULE_NAME Reco3D)

# Collect the headers
file(GLOB HEADERS "*.h")

# Remove LinkDef.h
list(FILTER HEADERS EXCLUDE REGEX ".*LinkDef.h$")

# However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cxx")

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# includes
include_directories(${LARCV_INCLUDE_DIR} $ENV{UBLARCVAPP_BASEDIR} ${LARLITE_INC_DIRS} ${OpenCV_INCLUDE_DIRS})

# Generate the dictionary
ROOT_GENERATE_DICTIONARY(G__${LIBNAME} ${HEADERS} LINKDEF LinkDef.h)

# Generate the shared library from the sources
add_library( ${LIBNAME} SHARED ${SOURCES} G__${LIBNAME}.cxx)
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_include_directories(${LIBNAME} PUBLIC ${LARCV_INCLUDE_DIR} ${LARLITE_INC_DIRS} ${OpenCV_INCLUDE_DIRS})
target_link_libraries(${LIBNAME} LArCVCoreBase LArCVCoreDataFormat opencv_core)

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
