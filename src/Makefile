CXX      :=  c++
MPICXX   :=  mpic++
CXXFLAGS :=  -O2 -Wall -Wchar-subscripts -g 

SOURCE   :=  main.o
SOURCE   +=  Constructors.o

SOURCE   += SpeedFunction.o 

SOURCE   += Integral.o section.o

SOURCE   += Transfer2.o 
SOURCE   += Transfer1.o 

SOURCE   += DataExchanger.o 

SOURCE   += MeshMpi.o 
SOURCE   += mesh/unstruct/UnstructMesh.o 

SOURCE   += mesh/unstruct/Polygon.o 
SOURCE   += mesh/unstruct/Prism.o mesh/unstruct/Tetrahedron.o mesh/unstruct/Hexahedron.o 

SOURCE   += mesh/unstruct/PhysicalFacet.o
SOURCE   += mesh/unstruct/MirrorFacet.o mesh/unstruct/GateFacet.o 
SOURCE   += mesh/unstruct/WallMaxwellFacet.o mesh/unstruct/MaxwellFacet.o 

SOURCE   += Printer.o 

SOURCE   +=  property_tree/property_tree.o
SOURCE   +=  base64/base64.o
SOURCE   +=  logger/logger.o

INCLUDES := -I .
INCLUDES += -I ../JsonCpp/include

LIBS     := -L ../JsonCpp/lib/ -ljsoncpp

.PHONY : all debug clean

all : kes

debug : CXXFLAGS += -DDEBUG -O0
debug : kes

clean :
	find . -iname "*.o" -exec rm '{}' ';'
	rm kes

kes : $(SOURCE)
	$(MPICXX) -o kes $(SOURCE) $(CXXFLAGS) $(LIBS)

MeshMpi.o : MeshMpi.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

DataExchanger.o : DataExchanger.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

base64/base64.o: base64/base64.cpp
	$(CXX) $(CXXFLAGS) -Wno-char-subscripts $(INCLUDES) -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@



