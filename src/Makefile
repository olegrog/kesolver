CC = mpic++
CCFLAGS = -O3 -I $(INCLUDE) -Wall -g
LIBXMLFLAGS = `xml2-config --cflags --libs`

PROGRAM_DIR = ./
SOURCE = main.o Constructors.o DataExchanger.o GateFacet.o Hexahedron.o Integral.o Loader.o MirrorFacet.o PhysicalFacet.o Polygon.o Printer.o Prism.o SpeedFunction.o Tetrahedron.o Transfer2.o Transfer1.o WallMaxwallFacet.o MaxwellFacet.o section.o

INCLUDE =  $(PROGRAM_DIR)

all : $(SOURCE)
	$(CC) -o kes $(SOURCE) $(CCFLAGS) $(LIBXMLFLAGS)
main.o : main.cpp
	$(CC) -c main.cpp $(CCFLAGS) $(LIBXMLFLAGS)
Loader.o : Loader.cpp
	$(CC) -c Loader.cpp $(CCFLAGS) $(LIBXMLFLAGS)
Constructors.o : Constructors.cpp ximesh.hpp ximesh_mixture.hpp ximesh_mix.hpp ximesh_rot.hpp ci.hpp ci_simple.hpp ci_mixture.hpp ci_mix.hpp maxwell.hpp macro.hpp stencil.hpp
	$(CC) -c Constructors.cpp $(CCFLAGS) $(LIBXMLFLAGS)
DataExchanger.o : DataExchanger.cpp
	$(CC) -c DataExchanger.cpp $(CCFLAGS) $(LIBXMLFLAGS)
Printer.o : Printer.cpp
	$(CC) -c Printer.cpp $(CCFLAGS) $(LIBXMLFLAGS)

Integral.o : Integral.cpp ci.hpp
	$(CC) -c Integral.cpp $(CCFLAGS)
section.o : section.cpp section.hpp
	$(CC) -c section.cpp $(CCFLAGS)

PhysicalFacet.o : PhysicalFacet.cpp
	$(CC) -c PhysicalFacet.cpp $(CCFLAGS)
GateFacet.o : GateFacet.cpp
	$(CC) -c GateFacet.cpp $(CCFLAGS)
MirrorFacet.o : MirrorFacet.cpp
	$(CC) -c MirrorFacet.cpp $(CCFLAGS)
WallMaxwallFacet.o : WallMaxwallFacet.cpp
	$(CC) -c WallMaxwallFacet.cpp $(CCFLAGS)
MaxwellFacet.o : MaxwellFacet.cpp
	$(CC) -c MaxwellFacet.cpp $(CCFLAGS)

Polygon.o : Polygon.cpp
	$(CC) -c Polygon.cpp $(CCFLAGS)
Tetrahedron.o : Tetrahedron.cpp
	$(CC) -c Tetrahedron.cpp $(CCFLAGS)
Prism.o : Prism.cpp
	$(CC) -c Prism.cpp $(CCFLAGS)
Hexahedron.o : Hexahedron.cpp
	$(CC) -c Hexahedron.cpp $(CCFLAGS)

SpeedFunction.o : SpeedFunction.cpp
	$(CC) -c SpeedFunction.cpp $(CCFLAGS)
XiMesh.o : XiMesh.cpp
	$(CC) -c XiMesh.cpp $(CCFLAGS)

Transfer1.o : Transfer1.cpp
	$(CC) -c Transfer1.cpp $(CCFLAGS) $(LIBXMLFLAGS)
Transfer2.o : Transfer2.cpp
	$(CC) -c Transfer2.cpp $(CCFLAGS) $(LIBXMLFLAGS)

