extern "C" {
#define ADD_
    // blas:
    #include <cblas_f77.h>

    // lapack:
    void dpotrf_( char *uplo, int *n, double *A, int *lda, int *info );
    void dtrtrs_( char *uplo, char *trans, char *diag, int *n, int *nrhs, double *A, int *lda,
         double *B, int *ldb, int *info );
    void dtrsm_( char *side, char *uplo, char *transA, char *diag, const int *m, const int *n, 
                 const double *alpha, const double *A, const int *lda, double *B, const int *ldb );
}

char boolToChar( bool value ) {
    return value ? 't' : 'n';
}

// double, general matrix multiply
void dgemm( bool transa, bool transb, int m, int n, int k, double alpha, double *A,
    int lda, double *B, int ldb, double beta, double *C, int ldc ) {
    char transachar = boolToChar( transa );
    char transbchar = boolToChar( transb );
    dgemm_(&transachar,&transbchar,&m,&n,&k, &alpha, A, &lda, B, &ldb, &beta, C, &ldc );    
}

// double, triangular, solve matrix
void dtrsm( bool XA, bool isUpper, bool transA, bool isUnitTriangular, int m, int n,
       double alpha, double *A, int lda, double *B, int ldb ) {
    char sideChar = XA ? 'R' : 'L';
    char isUpperChar = isUpper ? 'U' : 'L';
    char transAChar = transA ? 'T' : 'N';
    char isUnitTriangularChar = isUnitTriangular ? 'U' : 'N';
    dtrsm_( &sideChar, &isUpperChar, &transAChar, &isUnitTriangularChar, &m, &n,
        &alpha, A, &lda, B, &ldb );
}

// double, symmetric positive definite, triangular factorization (=cholesky)
int dpotrf( bool isUpper, int N, double *A, int lda ) {
    int info;
    char uplo = isUpper ? 'U' : 'L';
    dpotrf_( &uplo, &N, A, &lda, &info );
    return info;
}

// double, triangular, triangular solve
int dtrtrs( bool isUpper, bool transA, bool isUnitTriangular, int n, int nrhs, double *A, int lda,
             double *B, int ldb ) {
    int info;
    char isUpperChar = isUpper ? 'U' : 'L';
    char transChar = transA ? 'T' : 'N';
    char isUnitTriangularChar = isUnitTriangular ? 'U' : 'N';
    dtrtrs_( &isUpperChar, &transChar, &isUnitTriangularChar, &n, &nrhs, A, &lda, B, &ldb, &info );
    return info;
}
