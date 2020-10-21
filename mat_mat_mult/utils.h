#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double* generate_matrix(int n_row, int n_col)
{
  int i, j, idx;
  double* matrix = (double *)malloc(n_row*n_col*sizeof(double));

  for (i = 0; i < n_row; i++)
  {
    for (j = 0; j < n_col; j++)
    {
      idx = i*n_col + j;
      matrix[idx] = drand48();
      //matrix[idx] = (double)idx;
    }
  }

  return matrix;
}

double* reference_matrix(int n, double* A, double* B)
{
  int i, j, k;
  double* C = (double *)malloc(n*n*sizeof(double));

  for (k = 0; k < n; k++)
  {
    for (i = 0; i < n; i++)
    {
      double x = A[i*n+k];
      for (j = 0; j < n; j++)
	      C[i*n+j] += x * B[k*n+j];
    }
  }

  return C;
}

inline int fequal(double a, double b)
{
  double epsilon = 1e-10;
  if (fabs(a-b) < epsilon)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

inline void print_matrix(double* start, int n_row, int n_col)
{
  int i, j, idx;
  for (i = 0; i < n_row; i++)
  {
    printf("|  ");
    for (j = 0; j < n_col; j++)
    {
      idx = i*n_col + j;
      printf("%.0f  ", start[idx]);
    }
    printf("|\n");
  }
}

inline int check_matrix(int n, double* A, double* B)
{
  int i, j, idx;
  int error_counter=0;
  for (i = 0; i < n; i++)
  {
    for (j = 0; j < n; j++)
    {
      idx = i*n + j;
      if ( fequal(A[idx], B[idx]) == 0)
      {
        error_counter++;
        //printf("ERROR:A[%d]=%.16f\tB[%d]=%.16f\n", idx, A[idx], idx, B[idx]);
      }
    }
  }
  return error_counter;
}

inline void mat_mul(int n, double* A, double* B, double* C)
{
  int i, j, k;
  for (k = 0; k < n; k++)
  {
    for (i = 0; i < n; i++)
    {
      double x = A[i*n+k];
      for (j = 0; j < n; j++)
	      C[i*n+j] += x * B[k*n+j];
    }
  }
}

inline void mat_add(int n, double* A, double* B, double* C)
{
  int i, j;
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	      C[i*n+j] = A[i*n+j] + B[i*n+j];
    }
}

inline void mat_sub(int n, double* A, double* B, double* C)
{
  int i, j;
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
	      C[i*n+j] = A[i*n+j] - B[i*n+j];
    }
}
