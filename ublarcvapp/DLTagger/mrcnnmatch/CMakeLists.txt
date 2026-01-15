
# submodule name
set(SUBMODULE_NAME MRCNNMatch)

# this is a submodule of LARCVAPP_DLTagger
# for the source below, good practice is to get full paths, so we need
# the directory of target property
get_target_property(SOURCE_DIR LArCVApp_DLTagger SOURCE_DIR)

set(SOURCES
  ${SOURCE_DIR}/mrcnnmatch/MRCNNMatch.cxx
  ${SOURCE_DIR}/mrcnnmatch/MRCNNMatchConfig.cxx  
  ${SOURCE_DIR}/mrcnnmatch/MRCNNMatchTypes.cxx
  ${SOURCE_DIR}/mrcnnmatch/CropMaskCombo.cxx
  ${SOURCE_DIR}/mrcnnmatch/FeaturesMaskCombo.cxx
  ${SOURCE_DIR}/mrcnnmatch/Gen3DEndpoints.cxx
  ${SOURCE_DIR}/mrcnnmatch/GenGraphPoints.cxx  
  ${SOURCE_DIR}/mrcnnmatch/AStarMaskCombo.cxx
  )

# ADD THE SOURCES
target_sources( ${LIBNAME} PRIVATE ${SOURCES} )

# append the list of public headers for the target
set_property( TARGET ${LIBNAME}
  APPEND PROPERTY PUBLIC_HEADER
  ${SOURCE_DIR}/mrcnnmatch/MRCNNMatch.h
  ${SOURCE_DIR}/mrcnnmatch/MRCNNMatchConfig.h  
  ${SOURCE_DIR}/mrcnnmatch/MRCNNMatchTypes.h
  ${SOURCE_DIR}/mrcnnmatch/CropMaskCombo.h
  ${SOURCE_DIR}/mrcnnmatch/FeaturesMaskCombo.h
  ${SOURCE_DIR}/mrcnnmatch/Gen3DEndpoints.h
  ${SOURCE_DIR}/mrcnnmatch/GenGraphPoints.h  
  ${SOURCE_DIR}/mrcnnmatch/AStarMaskCombo.h
  )
  
