CC=g++
CFLAGS=-c -Wall -std=c++11
LDFLAGS=  -lblas 
SOURCES= $(wildcard *.cpp) 
OBJECTS=$(SOURCES:.cpp=.o)  
INCLUDES= -I. 

EXECUTABLE= gemm

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o


