

if (HAS_LARLITE)
  add_subdirectory(LArliteHandler)
  add_subdirectory(UBImageMod)
  add_subdirectory(UBWireTool)
  add_subdirectory(ubdllee)
  add_subdirectory(Filter)
endif()


if (LARCV_USE_GEO2D)
  add_subdirectory(LArOpenCVHandle)
endif()  

add_subdirectory(Stitchers)

