#add_subdirectory(LLCVProcessor)
add_subdirectory(dbscan)
#add_subdirectory(Reco3D)
#add_subdirectory(ContourTools)

if (HAS_LARLITE)
  add_subdirectory(LArliteHandler)
  add_subdirectory(UBImageMod)
  add_subdirectory(UBWireTool)
  add_subdirectory(ubdllee)
  add_subdirectory(Filter)
  add_subdirectory(MCTools)
  add_subdirectory(RecoTools)  
endif()


#if (LARCV_USE_GEO2D)
#  add_subdirectory(LArOpenCVHandle)
#endif()

#if (Boost_FOUND)
#  add_subdirectory(DLTagger)
#endif()
