

if (HAS_LARLITE)
  add_subdirectory(LArliteHandler)
  add_subdirectory(UBImageMod)
  add_subdirectory(ubdllee)
endif()


if (LARCV_USE_GEO2D)
  add_subdirectory(LArOpenCVHandle)
endif()  


