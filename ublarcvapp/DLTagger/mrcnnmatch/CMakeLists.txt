
# submodule name
set(SUBMODULE_NAME MRCNNMatch)

# However, the file(GLOB...) allows for wildcard additions:
set(SOURCES
  mrcnnmatch/MRCNNMatch.cxx
  mrcnnmatch/MRCNNMatchConfig.cxx  
  mrcnnmatch/MRCNNMatchTypes.cxx
  mrcnnmatch/CropMaskCombo.cxx
  mrcnnmatch/FeaturesMaskCombo.cxx
  mrcnnmatch/Gen3DEndpoints.cxx
  mrcnnmatch/GenGraphPoints.cxx  
  mrcnnmatch/AStarMaskCombo.cxx
  )

# ADD THE SOURCES
target_sources( ${LIBNAME} PRIVATE ${SOURCES} )

# append the list of public headers for the target
set_property( TARGET ${LIBNAME}
  APPEND PROPERTY PUBLIC_HEADER
  mrcnnmatch/MRCNNMatch.h
  mrcnnmatch/MRCNNMatchConfig.h  
  mrcnnmatch/MRCNNMatchTypes.h
  mrcnnmatch/CropMaskCombo.h
  mrcnnmatch/FeaturesMaskCombo.h
  mrcnnmatch/Gen3DEndpoints.h
  mrcnnmatch/GenGraphPoints.h  
  mrcnnmatch/AStarMaskCombo.h
  )
  
