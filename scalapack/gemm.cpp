#include <iostream>
#include <stdexcept>
#include <cstring>
#include <cmath>
#include <cblas.h>              /* Basic Linear Algebra I/O */
#include <chrono>
#include <limits.h>
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
    FILE * result;
    int n;
    int numthreads;
    Args( argc, argv ).arg("N", &n ).arg("numthreads", &numthreads ).go();
    openblas_set_num_threads( numthreads );

    int nprows = getRootFactor(P);
    int npcols = P / nprows;
//    if( p == 0 ) cout << "grid: " << nprows << " x " << npcols << endl;

    int system = blacs_get( -1, 0 );
    int grid = blacs_gridinit( system, true, nprows, npcols );
//    if( p == 0 ) cout << "system context " << system << " grid context: " << grid << endl;

    int myrow, mycol;
    blacs_gridinfo( grid, nprows, npcols, &myrow, &mycol );

    if( myrow >= nprows || mycol >= npcols ) {
        blacs_gridexit( grid );
        blacs_exit(0);
        exit(0);
    }

    // A     B       C
    // m x k k x n = m x n
    double * A;
    double * B;
    double * C; 
    double * C_ref; 

    if (p ==0){
      A = (double *)malloc(n * n * sizeof(double));
      B = (double *)malloc(n * n * sizeof(double));
      C = (double *)malloc(n * n * sizeof(double));
      C_ref = (double *)malloc(n * n * sizeof(double));
      result = fopen("result.txt", "a"); 
    }
    // nprows: process grid, number rows
    // npcols: process grid, number cols
    // myrow: process grid, our row
    // mycol: process grid, our col
    int m = n;
    int k = n;
    int nb=min(n/P,128);

    int mp = numroc( m, nb, myrow, 0, nprows ); // mp number rows A owned by this process
    int kp = numroc( k, nb, myrow, 0, nprows ); // kp number rows B owned by this process
    int kq = numroc( k, nb, mycol, 0, npcols ); // kq number cols A owned by this process
    int nq = numroc( n, nb, mycol, 0, npcols ); // nq number cols B owned by this process

    struct DESC desca, descb, descc;
    descinit( (&desca), m, k, nb, nb, 0, 0, grid, max(1, mp) );
    descinit( (&descb), k, n, nb, nb, 0, 0, grid, max(1, kp) );
    descinit( (&descc), m, n, nb, nb, 0, 0, grid, max(1, mp) );
    double *ipa = new double[desca.LLD_ * kq];
    double *ipb = new double[descb.LLD_ * nq];
    double *ipc = new double[descc.LLD_ * nq];


//    if( p == 0 ) cout << "created matrices" << endl;
    if( p == 0 ) {
      for(int i = 0; i < n * n; i++){
	A[i] = double(rand())/INT_MAX;
	B[i] = double(rand())/INT_MAX;
      }
    }
    double *work = new double[nb];
    if( n <=5 ) {
        pdlaprnt( n, n, ipa, 1, 1, &desca, 0, 0, "A", 6, work );
        pdlaprnt( n, n, ipb, 1, 1, &descb, 0, 0, "B", 6, work );
    }
    
    /// Scatter
    auto start_time =  std::chrono::steady_clock::now(); 
    int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
    for (int r = 0; r < n; r += nb, sendr=(sendr+1)%nprows) {
      sendc = 0;
      // Number of rows to be sent
      // Is this the last row block?
      int nr = nb;
      if (n-r < nb)
	  nr = n-r;
   
      for (int c = 0; c < n; c += nb, sendc=(sendc+1)%npcols) {
	  // Number of cols to be sent
	  // Is this the last col block?
	  int nc = nb;
	  if (n-c < nb)
	      nc = n-c;
   
	  if (p == 0) {
	      // Send a nr-by-nc submatrix to process (sendr, sendc)
	      Cdgesd2d(grid, nr, nc, A+n*c+r, n, sendr, sendc);
	      Cdgesd2d(grid, nr, nc, B+n*c+r, n, sendr, sendc);
	  }
   
	  if (myrow == sendr && mycol == sendc) {
	      // Receive the same data
	      // The leading dimension of the local matrix is nrows!
	      Cdgerv2d(grid, nr, nc, ipa+mp*recvc+recvr, mp, 0, 0);
	      Cdgerv2d(grid, nr, nc, ipb+mp*recvc+recvr, mp, 0, 0);
	      recvc = (recvc+nc)%kq;
	  }
   
      }
   
      if (myrow == sendr)
	  recvr = (recvr+nr)%mp;
    } 

//    Cpdgemr2d(n, n, Aseq, 0, 0, &descA_1x1, Apar, 0, 0, &descA_PxQ, grid); 

    MPI_Barrier( MPI_COMM_WORLD );

    pdgemm( false, false, m, n, k, 1,
                      ipa, 1, 1, &desca, ipb, 1, 1, &descb,
                      1, ipc, 1, 1, &descc );
    
    MPI_Barrier( MPI_COMM_WORLD );
    sendr = 0;
    for (int r = 0; r < n; r += nb, sendr=(sendr+1)%nprows) {
        sendc = 0;
        // Number of rows to be sent
        // Is this the last row block?
        int nr = nb;
        if (n-r < nb)
            nr = n-r;
 
        for (int c = 0; c < n; c += nb, sendc=(sendc+1)%npcols) {
            // Number of cols to be sent
            // Is this the last col block?
            int nc = nb;
            if (n-c < nb)
                nc = n-c;
 
            if (myrow == sendr && mycol == sendc) {
                // Send a nr-by-nc submatrix to process (sendr, sendc)
                Cdgesd2d(grid, nr, nc, ipc+mp*recvc+recvr, mp, 0, 0);
                recvc = (recvc+nc)%nq;
            }
 
            if (p == 0) {
                // Receive the same data
                // The leading dimension of the local matrix is nrows!
                Cdgerv2d(grid, nr, nc, C+n*c+r, n, sendr, sendc);
            }
 
        }
 
        if (myrow == sendr)
            recvr = (recvr+nr)%mp;
    }
    MPI_Barrier( MPI_COMM_WORLD );

    if (p ==0){
      auto end_time =  std::chrono::steady_clock::now(); 
      fprintf(result, "%d %d ", n, numthreads);
      fprintf(result, "%f ", std::chrono::duration_cast<std::chrono::milliseconds> (end_time - start_time).count() / 1000.0);
      
      start_time =  std::chrono::steady_clock::now(); 
      cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n,    n,    n, 1.0,   A,   n, B, n,  0.0, C_ref,  n);
      end_time =  std::chrono::steady_clock::now(); 
      fprintf(result, "%f ", std::chrono::duration_cast<std::chrono::milliseconds> (end_time - start_time).count() / 1000.0);
/*      for (int i =0; i < n; i++) {
	for (int j =0; j < n; j++)
	    printf("%f ", A[i * n + j]);
	printf("\n");
      }
      printf("\n");
      for (int i =0; i < n; i++) {
	for (int j =0; j < n; j++)
	    printf("%f ", B[i * n + j]);
	printf("\n");
      }
      printf("\n");
      for (int i =0; i < n; i++) {
	for (int j =0; j < n; j++)
	    printf("%f ", C[i * n + j]);
	printf("\n");
      }
      printf("\n");
*/      
      int error = 0;
      int count = 0;
      for (int i =0; i < n * n; i++){
	if (fabs(C[i] - C_ref[i]) > 0.1){
	  error = 1;
	  if (count <10){
	    printf("%f\t%f\n", C[i], C_ref[i]);
	    count++;
	  }
	}
      }
      if (error ==1)
	fprintf(result, "Failed\n");
      else
	fprintf(result, "Passed\n");
    }

    if (p ==0){
      free(A);
      free(B);
      free(C);
      free(C_ref);
      fclose(result);
    }
    free(ipa);
    free(ipb);
    free(ipc);

    blacs_gridexit( grid );
    blacs_exit(0);

    return 0;
}
