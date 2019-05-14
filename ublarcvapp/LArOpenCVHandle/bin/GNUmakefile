CC=g++
CFLAGS=-g
INCFLAGS =-I.
INCFLAGS+=-I$(shell ublarcvapp-config --incdir)
INCFLAGS+=-I$(shell larcv-config --incdir)
INCFLAGS+=$(shell laropencv-config --includes)
INCFLAGS+=$(shell geo2d-config --includes)
INCFLAGS+=-I$(shell larlite-config --incdir)
INCFLAGS+=-I$(shell larlite-config --incdir)/../UserDev
INCFLAGS+=-I$(OPENCV_INCDIR)
INCFLAGS+=$(shell root-config --cflags)


LDFLAGS=-L.
LDFLAGS+=$(shell ublarcvapp-config --libs) -lLArCVApp_UBImageMod -lLArCVApp_Filter
LDFLAGS+=$(shell larcv-config --libs) -lLArCVAppImageMod -lLArCVCoreCPPUtil
LDFLAGS+=$(shell laropencv-config --libs)
LDFLAGS+=$(shell geo2d-config --libs)
LDFLAGS+=$(shell larlite-config --libs)
LDFLAGS+=$(shell root-config --libs)


SOURCES=run.cxx
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=run

all: $(OBJECTS) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS)

$(OBJECTS): $(SOURCES)
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@

clean:
	rm *.o run
