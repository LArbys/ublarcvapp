set(MODULE_NAME Tagger)

# Collect the headers
file(GLOB HEADERS "*.h")

# Remove LinkDef.h
list(FILTER HEADERS EXCLUDE REGEX ".*LinkDef.h$")
 
# However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "*.cxx")

# library name
set(LIBNAME LArCVApp_${MODULE_NAME})
#message("laropencv libs: ${LAROPENCV_LIBS}")

# larlite includes
set(LARLITE_INCLUDES $ENV{LARLITE_USERDEVDIR} $ENV{LARLITE_COREDIR} $ENV{LARLITE_USERDEVDIR}/BasicTool )
#message("larlite inc: ${LARLITE_INCLUDES}")

# larlite libs used
#set(LARLITE_LIBS_USED SelectionTool_LEEPreCuts

# larcv libs used
set(LARCV_LIBS_USED LArCVCoreBase LArCVCoreDataFormat LArCVCoreCVUtil)

# ublarcvapp libraries
set(UBLARCVAPP_LIBS_USED LArCVApp_ContourTools LArCVApp_UBImageMod)

# includes
include_directories(${LARCV_INCLUDE_DIR} $ENV{UBLARCVAPP_BASEDIR} $ENV{LAROPENCV_BASEDIR} ${OpenCV_INCLUDE_DIRS} ${GEO2D_INC_DIR} ${LARLITE_INCLUDES} )

# Generate the dictionary
ROOT_GENERATE_DICTIONARY(G__${LIBNAME} ${HEADERS} LINKDEF LinkDef.h)

# Generate the shared library from the sources
add_library( ${LIBNAME} SHARED ${SOURCES} G__${LIBNAME}.cxx)
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_include_directories(${LIBNAME} PUBLIC ${LARCV_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS} ${LARLITE_INCLUDES} )
target_link_libraries(${LIBNAME} ${LARCV_LIBS_USED} ${UBLARCVAPP_LIBS_USED} ${LAROPENCV_LIBS} ${GEO2D_LIBS} opencv_core)

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
