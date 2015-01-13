export OPENBLAS_NUM_THREADS=2
CC=mpic++.mpich2
CFLAGS=-c -Wall -std=c++11 
LDFLAGS= -L../../scalapackInstall/scalapack-2.0.2/ -lscalapack -L/opt/OpenBLAS/lib -lopenblas -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.8/ -lgfortran
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


