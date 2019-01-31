set(MODULE_NAME LArOpenCVHandle)

# Collect the headers
file(GLOB HEADERS "*.h")

# Remove LinkDef.h
list(FILTER HEADERS EXCLUDE REGEX ".*LinkDef.h$")
 
# However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cxx")

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})

# includes
include_directories(${LARCV_INCLUDE_DIR} ${GEO2D_INC_DIR} ${LAROPENCV_INC_DIR} ${LARLITE_INC_DIRS})

# Generate the dictionary
ROOT_GENERATE_DICTIONARY(${LIBNAME}Dict ${HEADERS} LINKDEF LinkDef.h)

# Generate the shared library from the sources
add_library( ${LIBNAME} SHARED ${SOURCES} ${LIBNAME}Dict.cxx)
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_include_directories(${LIBNAME} PUBLIC ${LARCV_INCLUDE_DIR} ${GEO2D_INC_DIR} ${LARLITE_INC_DIRS})
target_link_libraries(${LIBNAME}
  LArCVCoreBase LArCVCoreDataFormat LArCVCoreProcessor 
  ${LAROPENCV_LIBS}
  ${GEO2D_LIBS} )

# install libraries and headers
install(TARGETS ${LIBNAME}
  EXPORT ublarcvapp
  LIBRARY DESTINATION lib
  PUBLIC_HEADER DESTINATION include/ublarcvapp/${MODULE_NAME})

# install files needed for dictionary use in CINT/CLING/PYTHON
install( FILES
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}Dict_rdict.pcm
  ${CMAKE_CURRENT_BINARY_DIR}/lib${LIBNAME}Dict.rootmap
  DESTINATION lib  )
