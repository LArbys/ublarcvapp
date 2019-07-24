CXX=g++
CFLAGS=-g -fPIC
INCLUDES =  `root-config --cflags`
INCLUDES += `larcv-config --includes`
INCLUDES += `larlite-config --includes` -I${LARLITE_USERDEVDIR} -I${LARLITE_USERDEVDIR}/BasicTool
INCLUDES += `ublarcvapp-config --includes`-I${UBLARCVAPP_INCDIR}/ublarcvapp/Tagger

LDLIBS =
LDLIBS += `ublarcvapp-config --libs`
LDLIBS += -lLArCVApp_UBImageMod -lLArCVApp_Tagger -lLArCVApp_LArliteHandler -lLArCVApp_UBWireTool
LDLIBS += -lLArCVApp_dbscan
LDLIBS += `larcv-config --libs`  -lLArCVCoreCVUtil
LDLIBS += `laropencv-config --libs`
LDLIBS += `geo2d-config --libs`
LDLIBS += `larlite-config --libs` -lSelectionTool_LEEPreCuts -lBasicTool_FhiclLite
LDLIBS +=  `root-config --libs`
LDLIBS += -lopencv_imgcodecs
LDLIBS += -lstdc++

all: run_ubdl_tagger

%: %.cxx
	$(CXX) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LDLIBS)

clean:
	rm run_ubdl_tagger
