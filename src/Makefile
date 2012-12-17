CXX       =  c++
MPICXX    =  mpic++
CXXFLAGS  =  -O3 -Wall -Wchar-subscripts -g -msse2

SOURCE    =  main.o
SOURCE   +=  Constructors.o
SOURCE   +=  property_tree/property_tree.o
SOURCE   +=  base64/base64.o

INCLUDES := -I ../JsonCpp/include
LIBS     := -L ../JsonCpp/lib/ -ljsoncpp

all : $(SOURCE)
	$(MPICXX) -o kes $(SOURCE) $(CXXFLAGS) $(LIBS)

main.o : main.cpp
	$(MPICXX) -c $< -o $@ $(CXXFLAGS) $(INCLUDES)

base64/base64.o: base64/base64.cpp
	$(CXX) $(CXXFLAGS) -Wno-char-subscripts $(INCLUDES) -c $< -o $@

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

