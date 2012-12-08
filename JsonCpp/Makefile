CXX      := g++
CXXFLAGS := -c -Wall -O2
AR       := ar 
ARFLAGS  := rcs 

SOURCES  := src/json_value.cpp
SOURCES  += src/json_reader.cpp
SOURCES  += src/json_writer.cpp

INCLUDES := -I include/
LIBS     := 

TARGET   := lib/libjsoncpp.a

OBJECTS=$(SOURCES:.cpp=.o)

all: $(TARGET) $(OBJECTS) 
	
$(TARGET): $(OBJECTS) 
	$(AR) $(ARFLAGS) $@ $(OBJECTS) 

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) $< -o $@

clean:
	rm $(OBJECTS) $(TARGET)
