/**
 * svd_demo.cpp
 * 
 * Given that you put version 3.5.0 into /opt/lapack/ compile this with: 
 * g++ svd_demo.cpp -I"/opt/lapack/lapack-3.5.0/lapacke/include" \
 *   -L"/opt/lapack/lapack-3.5.0" -llapacke -llapack -lblas -lcblas
 * The order of included libraries is important!
 */

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <cblas.h>
#include <lapacke.h>

using namespace std;

typedef double value;

/** Column major style! */
string matrix2string(int m, int n, value* A)
{
  ostringstream oss;
  for (int j=0;j<m;j++)
  {
    for (int k=0;k<n;k++)
    {
      oss << A[j+k*m] << "\t";
    }
    oss << endl;
  }
  return oss.str();
}

int main(int argc, char** argv)
{
  //> Part 1. Decomposition. -----------------------------------------
  char jobu  = 'A'; // Return the complete matrix U
  char jobvt = 'A'; // Return the complete matrix VT
  int mA = 2;
  int nA = 3;
  int lda = 2;
  int ldu = 2;
  int ldvt = 3;
  int lwork = 81;
  int info = 0;
  value* A = (value*)malloc(mA*nA*sizeof(value));
  value* U = (value*)malloc(mA*mA*sizeof(value));
  value* VT = (value*)malloc(nA*nA*sizeof(value));
  value* Svec = (value*)malloc(3*sizeof(value));
  value* work = (value*)malloc(lwork*sizeof(value));

  A[0] = 1; A[2] = 2; A[4] = 4;
  A[1] = 0; A[3] = 0; A[5] = 4;

  cout << "Matrix A (will be overwritten, as is documented):" << endl <<
    matrix2string(mA,nA,A);

  // Citing lapacke.h
  lapack_int LAPACKE_dgesvd(int matrix_order, char jobu, char jobvt,
     lapack_int m, lapack_int n, double* a,
     lapack_int lda, double* s, double* u, lapack_int ldu,
     double* vt, lapack_int ldvt, double* superb);

  info = LAPACKE_dgesvd(LAPACK_COL_MAJOR, jobu, jobvt, mA, nA, A, lda, Svec, U, ldu, VT, ldvt, work);
  cout << "Ran dgesvd. Let's see ..." << endl <<
    "U:" << endl << matrix2string(mA,mA,U) <<
    "Svec:" << endl << matrix2string(1,nA,Svec) <<
    "VT:" << endl << matrix2string(nA,nA,VT) <<
    "Info Code: " << info << endl << endl <<
    "All is well." << endl;
  //< ----------------------------------------------------------------
  //> Part 2. Checking the result. -----------------------------------
  value* S = (value*)malloc(mA*nA*sizeof(value));
  S[0] = Svec[0]; S[2] = 0      ; S[4] = 0      ;
  S[1] = 0      ; S[3] = Svec[1]; S[5] = 0      ;

  // Citing cblas.h
  // void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
  //   const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
  //   const int K, const double alpha, const double *A,
  //   const int lda, const double *B, const int ldb,
  //   const double beta, double *C, const int ldc);

  // work := S*VT; (2x3)=(2x3)*(3x3)
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mA,nA,nA,1,S,lda,VT,ldvt,0,work,lda)    ;
  cout << "Step 1: S*VT" << endl << matrix2string(2,3,work);

  // A := U*work; (2x2)*(2x3)
  cblas_dgemm(CblasColMajor,CblasNoTrans,CblasNoTrans,mA,nA,mA,1,U,ldu,work,lda,0,A,lda);
  cout << "A := U*S*VT:" << endl << matrix2string(mA,nA,A) << endl;
  //< ----------------------------------------------------------------
  free(A); free(U); free(VT); free(Svec); free(work); free(S);
  return EXIT_SUCCESS;
}
