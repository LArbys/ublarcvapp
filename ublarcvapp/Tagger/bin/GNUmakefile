CXX=g++
CFLAGS=-g -fPIC
INCLUDES =  `root-config --cflags`
INCLUDES += `larcv-config --includes`
INCLUDES += `larlite-config --includes` -I${LARLITE_USERDEVDIR}
INCLUDES += `ublarcvapp-config --includes`

LDLIBS =
LDLIBS += `ublarcvapp-config --libs` -lLArCVApp_UBImageMod -lLArCVApp_Tagger -lLArCVApp_LArliteHandler
LDLIBS += `larcv-config --libs`  -lLArCVCoreCVUtil
LDLIBS += `laropencv-config --libs`
LDLIBS += `geo2d-config --libs`
LDLIBS += `larlite-config --libs`
LDLIBS +=  `root-config --libs`
LDLIBS += -lstdc++

all: run_ubdl_tagger

%: %.cxx
	$(CXX) $(CFLAGS) $(INCLUDES) $^ -o $@ $(LDLIBS)

clean:
	rm run_ubdl_tagger
