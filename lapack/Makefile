export OPENBLAS_NUM_THREADS=2
CC=g++
CFLAGS=-c -Wall -std=c++11 
LDFLAGS=  -L/opt/OpenBLAS/lib -lopenblas -llapack
SOURCES= $(wildcard *.cpp) 
OBJECTS=$(SOURCES:.cpp=.o)  
INCLUDES=  -I. -I/opt/OpenBLAS/include

EXECUTABLE= gemm

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) $(LDFLAGS) -o $@ 

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o


