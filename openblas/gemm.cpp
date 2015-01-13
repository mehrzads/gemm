
#include <stdio.h>              /* I/O lib         ISOC     */
#include <stdlib.h>             /* Standard Lib    ISOC     */
#include <limits.h>             /* Standard Lib    ISOC     */
#include <cblas.h>              /* Basic Linear Algebra I/O */
#include <chrono>
#define eps 0

void matrixMultiply(int N, double * A, double * B, double * C){
  for (int i = 0; i < N; i++)
    for (int j = 0; j < N ; j++){
      double sum = 0;
      for (int k = 0; k < N ; k++)
	sum += A[i *N + k] * B[k * N + j];
      C[i * N + j] = sum;
    }
  return;

}

int main(int argc, char **argv) {
  int N= 8192;
  double * A = (double *)malloc( N * N * sizeof(double));
  double * B = (double *)malloc( N * N * sizeof(double));
  double * C = (double *)malloc( N * N * sizeof(double));
  double * C_ref = (double *)malloc( N * N * sizeof(double));

  for (int i = 0; i < N * N; i ++){
    A[i] = double(rand())/INT_MAX;
    B[i] = double(rand())/INT_MAX;
    C[i] = 0;
    C_ref[i] = 0;
  }
  printf("%s\n", openblas_get_config());
  openblas_set_num_threads(2);
  auto start_time =  std::chrono::steady_clock::now(); 
            /* row_order      transform     transform     rowsA colsB K  alpha  a  lda  b  ldb beta c   ldc */
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N,    N,    N, 1.0,   A,   N, B, N,  0.0, C,  N);
  auto end_time =  std::chrono::steady_clock::now(); 

  printf("blas takes %f Seconds\n", std::chrono::duration_cast<std::chrono::milliseconds> (end_time - start_time).count() / 1000.0);

  start_time =  std::chrono::steady_clock::now(); 
  matrixMultiply(N, A, B, C_ref);
  end_time =  std::chrono::steady_clock::now(); 

  printf("naive takes %f Seconds\n", std::chrono::duration_cast<std::chrono::milliseconds> (end_time - start_time).count() / 1000.0);

  int ferror =0;
  int count = 0;
  for (int i = 0; i < N * N; i++){
    if (abs(C[i] - C_ref[i]) > eps)
    {
      ferror = 1;
      if (count < 10){
	printf("%d\t%f\t%f\n", i, C[i], C_ref[i]);
	count++;
      }
    }
  }
  if (ferror !=0)
    printf("Failed\n");
  else
    printf("Passed\n");


  free(A);
  free(B);
  free(C);
  free(C_ref);

}
