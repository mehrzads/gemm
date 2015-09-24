#!/bin/bash
#mpirun.mpich2 -n 1 -f host_file ./gemm 8192 8
mpirun.mpich2 -n 1  ./gemm 16384 8
#mpirun.mpich2 -n 1  ./gemm 16384 4


# $Id: $


