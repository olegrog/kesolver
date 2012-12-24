CXX      :=  c++
MPICXX   :=  mpic++ -DMPI_ON
CXXFLAGS :=  -O2 -Wall -Wchar-subscripts -g 

SOURCE   :=  main.c
SOURCE   +=  Constructors.c

SOURCE   += SpeedFunction.c 

SOURCE   += Integral.c section.c

SOURCE   += Transfer2.c 
SOURCE   += Transfer1.c 

SOURCE   += MeshSingle.c 
SOURCE   += mesh/unstruct/UnstructMesh.c 

SOURCE   += mesh/unstruct/Polygon.c 
SOURCE   += mesh/unstruct/Prism.c mesh/unstruct/Tetrahedron.c mesh/unstruct/Hexahedron.c 

SOURCE   += mesh/unstruct/PhysicalFacet.c
SOURCE   += mesh/unstruct/MirrorFacet.c mesh/unstruct/GateFacet.c 
SOURCE   += mesh/unstruct/WallMaxwellFacet.c mesh/unstruct/MaxwellFacet.c 

SOURCE   += Printer.c 

SOURCE   +=  property_tree/property_tree.c
SOURCE   +=  base64/base64.c
SOURCE   +=  logger/logger.c

MPISOURCE := MeshMpi.c 
MPISOURCE += DataExchanger.c 

OBJ    = $(SOURCE:.c=.o)
MPIOBJ = $(MPISOURCE:.c=.o)

INCLUDES := -I .
INCLUDES += -I ../JsonCpp/include

LIBS     := -L ../JsonCpp/lib/ -ljsoncpp

.PHONY : all debug clean

all : kes kes-mpi

debug : CXXFLAGS += -DDEBUG -O0
debug : kes

clean :
	find . -iname "*.o" -exec rm '{}' ';'
	rm kes kes-mpi

kes : $(OBJ)
	$(CXX) -o kes $(OBJ) $(CXXFLAGS) $(LIBS)

kes-mpi : $(OBJ) $(MPIOBJ)
	$(MPICXX) -o kes-mpi $(OBJ) $(MPIOBJ) $(CXXFLAGS) $(LIBS)

MeshMpi.o : MeshMpi.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

DataExchanger.o : DataExchanger.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

base64/base64.o: base64/base64.cpp
	$(CXX) $(CXXFLAGS) -Wno-char-subscripts $(INCLUDES) -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@



