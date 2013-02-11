CXX      :=  c++
MPICXX   :=  mpic++ -DMPI_ON
CXXFLAGS :=  -O2 -Wall -g -pg -funroll-loops

SOURCE   :=  main.cpp
SOURCE   +=  Constructors.cpp

SOURCE   += SpeedFunction.cpp

SOURCE   += Integral.cpp section.cpp

SOURCE   += Transfer2.cpp
SOURCE   += Transfer1.cpp

SOURCE   += MeshSingle.cpp
SOURCE   += mesh/unstruct/UnstructMesh.cpp

SOURCE   += mesh/unstruct/Polygon.cpp
SOURCE   += mesh/unstruct/Prism.cpp mesh/unstruct/Tetrahedron.cpp mesh/unstruct/Hexahedron.cpp

SOURCE   += mesh/unstruct/PhysicalFacet.cpp
SOURCE   += mesh/unstruct/MirrorFacet.cpp mesh/unstruct/GateFacet.cpp
SOURCE   += mesh/unstruct/WallMaxwellFacet.cpp mesh/unstruct/MaxwellFacet.cpp

SOURCE   += Printer.cpp

SOURCE   +=  property_tree/property_tree.cpp
SOURCE   +=  base64/base64.cpp
SOURCE   +=  logger/logger.cpp

MPISOURCE := MeshMpi.cpp
MPISOURCE += DataExchanger.cpp

OBJ    = $(SOURCE:.cpp=.o)
MPIOBJ = $(MPISOURCE:.cpp=.o)

INCLUDES := -I .
INCLUDES += -I ../JsonCpp/include

LIBS     := -L ../JsonCpp/lib/ -ljsoncpp

DEP    = $(SOURCE:.cpp=.d)
MPIDEP = $(MPISOURCE:.cpp=.d)

.PHONY : all debug clean

all : kes kes-mpi

debug : CXXFLAGS += -DDEBUG -O0
debug : kes

clean :
	find . -iname "*.[od]" -exec rm '{}' ';'
	rm kes kes-mpi

kes : $(OBJ)
	$(CXX) -o kes $(OBJ) $(CXXFLAGS) $(LIBS)

kes-mpi : $(OBJ) $(MPIOBJ)
	$(MPICXX) -o kes-mpi $(OBJ) $(MPIOBJ) $(CXXFLAGS) $(LIBS)

-include $(DEP)
-include $(MPIDEP)

%.d: %.cpp
	@set -e; rm -f $@; \
	$(CXX) -MM $(CXXFLAGS) $(INCLUDES) $< > $@.$$$$; \
	sed 's,\($(notdir $*)\)\.o[ :]*,$(dir $*)\1.o $@ : ,g' < $@.$$$$ > $@; \
	rm -f $@.$$$$

MeshMpi.o : MeshMpi.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

DataExchanger.o : DataExchanger.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

base64/base64.o: base64/base64.cpp
	$(CXX) $(CXXFLAGS) -Wno-char-subscripts $(INCLUDES) -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

