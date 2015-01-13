#include <iostream>
#include <stdexcept>
#include <cstring>
#include <cmath>
using namespace std;

#include "mpi.h"

#include "utils/NanoTimer.h"
#include "utils/stringhelper.h"
#include "args.h"
#include "scalapack.h"

extern "C" {
    void openblas_set_num_threads(int num_threads);
}

int getRootFactor( int n ) {
    for( int t = sqrt(n); t > 0; t-- ) {
        if( n % t == 0 ) {
            return t;
        }
    }
    return 1;
}

// conventions:
// M_ by N_ matrix block-partitioned into MB_ by NB_ blocks, then
// distributed according to 2d block-cyclic scheme

// based on http://acts.nersc.gov/scalapack/hands-on/exercise3/pspblasdriver.f.html

int main( int argc, char *argv[] ) {
    int p, P;
    blacs_pinfo( &p, &P );
//    mpi_print( toString(p) + " / " + toString(P) );

    int n;
    int numthreads;
    int its;
    int blocksize;
    Args( argc, argv ).arg("N", &n ).arg("num iterations", &its ).arg("numthreads", &numthreads ).arg("blocksize", &blocksize).go();
    openblas_set_num_threads( numthreads );

    int nprows = getRootFactor(P);
    int npcols = P / nprows;
    if( p == 0 ) cout << "grid: " << nprows << " x " << npcols << endl;

    int system = blacs_get( -1, 0 );
    int grid = blacs_gridinit( system, true, nprows, npcols );
    if( p == 0 ) cout << "system context " << system << " grid context: " << grid << endl;

    int myrow, mycol;
    blacs_gridinfo( grid, nprows, npcols, &myrow, &mycol );
//    mpi_print("grid, me: " + toString(myrow) + ", " + toString(mycol) );

    if( myrow >= nprows || mycol >= npcols ) {
//        mpi_print("not needed, exiting");
        blacs_gridexit( grid );
        blacs_exit(0);
        exit(0);
    }

    // A     B       C
    // m x k k x n = m x n
    // nb: blocksize

    // nprows: process grid, number rows
    // npcols: process grid, number cols
    // myrow: process grid, our row
    // mycol: process grid, our col
    int m = n;
    int k = n;
//    int nb = min(n,128); // nb is column block size for A, and row blocks size for B
    int nb=min(n/P,128);

    int mp = numroc( m, nb, myrow, 0, nprows ); // mp number rows A owned by this process
    int kp = numroc( k, nb, myrow, 0, nprows ); // kp number rows B owned by this process
    int kq = numroc( k, nb, mycol, 0, npcols ); // kq number cols A owned by this process
    int nq = numroc( n, nb, mycol, 0, npcols ); // nq number cols B owned by this process
//    mpi_print( "mp " + toString(mp) + " kp " + toString(kp) + " kq " + toString(kq) + " nq " + toString(nq) );

    struct DESC desca, descb, descc;
    descinit( (&desca), m, k, nb, nb, 0, 0, grid, max(1, mp) );
    descinit( (&descb), k, n, nb, nb, 0, 0, grid, max(1, kp) );
    descinit( (&descc), m, n, nb, nb, 0, 0, grid, max(1, mp) );
//    mpi_print( "desca.LLD_ " + toString(desca.LLD_) + " kq " + toString(kq) );
    double *ipa = new double[desca.LLD_ * kq];
    double *ipb = new double[descb.LLD_ * nq];
    double *ipc = new double[descc.LLD_ * nq];

    for( int i = 0; i < desca.LLD_ * kq; i++ ) {
        ipa[i] = p;
    }
    for( int i = 0; i < descb.LLD_ * nq; i++ ) {
        ipb[i] = p;
    }

    if( p == 0 ) cout << "created matrices" << endl;
    double *work = new double[nb];
    if( n <=5 ) {
        pdlaprnt( n, n, ipa, 1, 1, &desca, 0, 0, "A", 6, work );
        pdlaprnt( n, n, ipb, 1, 1, &descb, 0, 0, "B", 6, work );
    }

    NanoTimer timer;
    for( int it = 0; it < its; it++ ) {
        pdgemm( false, false, m, n, k, 1,
                      ipa, 1, 1, &desca, ipb, 1, 1, &descb,
                      1, ipc, 1, 1, &descc );
        MPI_Barrier( MPI_COMM_WORLD );
        if( p == 0 ) timer.toc("it " + toString(it) + " pdgemm");
    }

    blacs_gridexit( grid );
    blacs_exit(0);

    return 0;
}
