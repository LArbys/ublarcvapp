set(MODULE_NAME LLCVProcessor)

# DEFINE THE TARGET
set( LIBNAME LArCVApp_${MODULE_NAME} )

set( HEADERS LLCVProcessDriver.h
  LLCVProcessBase.h
  LLCVProcessTest.h
  LArCVProcessTest.h )

add_library( ${LIBNAME} SHARED
  LLCVProcessDriver.cxx
  LLCVProcessBase.cxx
  LLCVProcessTest.cxx
  LArCVProcessTest.cxx )

# create a property to hold our headers
set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")

# DEFINE SOME VARS THAT MIGHT BE USEFUL FOR SUBMODULES

# larlite includes
set(LARLITE_INCLUDES $ENV{LARLITE_BASEDIR} $ENV{LARLITE_USERDEVDIR} $ENV{LARLITE_COREDIR} $ENV{LARLITE_USERDEVDIR}/BasicTool )

# larcv libs used
set(LARCV_LIBS_USED LArCVCoreBase LArCVCoreDataFormat LArCVCoreCVUtil)

# ublarcvapp libraries
set(UBLARCVAPP_LIBS_USED LArCVApp_LArliteHandler)

#message("headers: ${HEADERS}")
#message("project source dir: ${PROJECT_SOURCE_DIR}")
#message("eigen3: ${EIGEN3_INCLUDE_DIRS}")

get_target_property(ALL_HEADERS ${LIBNAME} PUBLIC_HEADER)

target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  ${LARCV_INCLUDE_DIR}
  ${LARLITE_INCLUDES}
  PRIVATE
  ${PROJECT_SOURCE_DIR}
  )
target_link_libraries(${LIBNAME} ${LARLITE_LIBS_USED} ${LARCV_LIBS_USED} ${UBLARCVAPP_LIBS_USED})

include_directories( ${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/ublarcvapp/${MODULE_NAME} ${LARCV_INCLUDE_DIR} ${LARLITE_INCLUDES} )
ROOT_GENERATE_DICTIONARY(G__${LIBNAME} ${HEADERS} LINKDEF LinkDef.h)
target_sources( ${LIBNAME} PRIVATE G__${LIBNAME}.cxx )

#message("public headers: ${ALL_HEADERS}")

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
