MPICXX   :=  mpic++.openmpi
CXX      :=  g++
CXXFLAGS :=  -std=c++11 -O3 -g3 -Wall -march=native -mno-avx # -fprofile-use

SOURCE   :=  Constructors.cpp

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
SOURCE   += mesh/unstruct/AlphaFacet.cpp

SOURCE   += Printer.cpp

SOURCE   +=  property_tree/property_tree.cpp
SOURCE   +=  base64/base64.cpp
SOURCE   +=  logger/logger.cpp

MPISOURCE := MeshMpi.cpp
MPISOURCE += DataExchanger.cpp

MPIOBJ_ALL := $(SOURCE:.cpp=.o)

SOURCE    +=  main.cpp
MPISOURCE +=  main-mpi.cpp

OBJ        := $(SOURCE:.cpp=.o)
MPIOBJ     := $(MPISOURCE:.cpp=.o)
MPIOBJ_ALL += $(MPISOURCE:.cpp=.o)

INCLUDES := -I .
INCLUDES += -I ../JsonCpp/include

LIBS     := -L ../JsonCpp/lib/ -ljsoncpp

DEP    := $(SOURCE:.cpp=.d)
MPIDEP := $(MPISOURCE:.cpp=.d)

.PHONY : all debug clean

all : kes kes-mpi

debug : CXXFLAGS += -DDEBUG -O0
debug : kes

clean :
	find . -iname "*.[od]*" -exec rm -r {} \;
	rm -f kes kes-mpi

kes : $(OBJ)
	$(CXX) -o kes $(OBJ) $(CXXFLAGS) $(LIBS)

kes-mpi : $(MPIOBJ_ALL)
	$(MPICXX) -o kes-mpi $(MPIOBJ_ALL) $(CXXFLAGS) $(LIBS)

-include $(DEP)
-include $(MPIDEP)

define make_dep
	@set -e; rm -f $1; \
	$3 -MM $(CXXFLAGS) $(INCLUDES) $< > $1.$$$$; \
	sed 's,\($(notdir $2)\)\.o[ :]*,$(dir $2)\1.o $1 : ,g' < $1.$$$$ > $1; \
	rm -f $1.$$$$
endef

$(DEP): %.d: %.cpp
	$(call make_dep,$@,$*,$(CXX))

$(MPIDEP): %.d: %.cpp
	$(call make_dep,$@,$*,$(MPICXX))

$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(MPIOBJ): %.o: %.cpp
	$(MPICXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

