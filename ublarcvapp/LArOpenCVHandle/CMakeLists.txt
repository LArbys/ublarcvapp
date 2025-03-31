set(MODULE_NAME LArOpenCVHandle)

set(LIBNAME LArCVApp_${MODULE_NAME})

# Collect the headers
set( HEADERS
  BadRegionAna.h
  BlankImage.h
  CosmicPixelAna.h
  CosmicTrackAna.h
  CosmicVertexAna.h
  DeadWireAna.h
  GenDeadChannelMap.h
  GoodcROIFilter.h
  ImageROICheckMC.h
  InterWriteOutLC.h
  LArbysImage.h
  LArbysImageExtract.h
  LArbysImageMC.h
  LArbysImageMaker.h
  LArbysLoader.h
  LArbysUtils.h
  OperateImage.h
  PGraphTruthMatch.h
  PIDImageMaker.h
  PostTagger.h
  PreProcessor.h
  ROIAna.h
  SSNetEventAna.h
  SSNetTruthAna.h
  SegmentFilter.h
  SingleROIFaker.h
  VertexAna.h
  VertexFilter.h
  VertexInROI.h
)

add_library( ${LIBNAME} SHARED
  BadRegionAna.cxx
  BlankImage.cxx
  CosmicPixelAna.cxx
  CosmicTrackAna.cxx
  CosmicVertexAna.cxx
  DeadWireAna.cxx
  GenDeadChannelMap.cxx
  GoodcROIFilter.cxx
  ImageROICheckMC.cxx
  InterWriteOutLC.cxx
  LArbysImage.cxx
  LArbysImageExtract.cxx
  LArbysImageMC.cxx
  LArbysImageMaker.cxx
  LArbysLoader.cxx
  LArbysUtils.cxx
  OperateImage.cxx
  PGraphTruthMatch.cxx
  PIDImageMaker.cxx
  PostTagger.cxx
  PreProcessor.cxx
  ROIAna.cxx
  SSNetEventAna.cxx
  SSNetTruthAna.cxx
  SegmentFilter.cxx
  SingleROIFaker.cxx
  VertexAna.cxx
  VertexFilter.cxx
  VertexInROI.cxx
)  

set_target_properties(${LIBNAME} PROPERTIES PUBLIC_HEADER "${HEADERS}")
target_include_directories(${LIBNAME}
  PUBLIC
  $<INSTALL_INTERFACE:include/ublarcvapp/${MODULE_NAME}>
  PRIVATE
  ${PROJECT_SOURCE_DIR}
  ${OpenCV_INCLUDE_DIRS}  
  ${LARLITE_INCLUDE_DIR}
  ${GEO2D_INC_DIR}
  ${LAROPENCV_INC_DIR}
  ${LARCV_INCLUDE_DIR}  
)
target_link_libraries( ${LIBNAME} PUBLIC
  opencv_imgcodecs
  opencv_core  
  LArCVCoreBase
  LArCVCoreDataFormat
  LArCVCoreCVUtil
  ${GEO2D_LIBS}  
  ${LAROPENCV_LIBS}  
  LArCVApp_Reco3D  
)
include_directories( ${CMAKE_CURRENT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}
  ${OpenCV_INCLUDE_DIRS}  
  ${LARLITE_INCLUDE_DIR}
  ${GEO2D_INC_DIR}
  ${LAROPENCV_INC_DIR}
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
