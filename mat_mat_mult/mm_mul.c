/*********************************************************************
Matrix-matrix multiplication with Fox's algorithm.

*********************************************************************/

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"

/*Uncomment to compare against reference matrix*/
//#define REF

int main(int argc, char *argv[])
{
  /*MPI_variables-------------------------------------------------------------*/
  int world_id, row_id, col_id, nproc, i, j, k, n_row, n_col;
  int count, blocklen, stride, idx;
  int coords[2], pos[2], reorder = 1, ndim = 2, nodes_per_dim[2] = {0,0}, periods[2] = {0,0};
  MPI_Comm proc_grid, proc_row, proc_col;
  MPI_Datatype block;
  MPI_Request request;
  MPI_Status status;
  double timer, avg_time, max_time;
  double *A, *B, *C, *Asub, *Bsub, *Csub, *Cref, *Ainc, *Binc;
  n_row = atoi(argv[1]);
  n_col = n_row;
  /*--------------------------------------------------------------------------*/

  MPI_Init(&argc, &argv); /* Initialize MPI */

  /*Create a virtual 2D-grid topology-----------------------------------------*/
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  nodes_per_dim[0] = sqrt(nproc);
  nodes_per_dim[1] = sqrt(nproc);
  MPI_Dims_create(nproc, ndim, nodes_per_dim);
  MPI_Cart_create(MPI_COMM_WORLD, ndim, nodes_per_dim, periods, reorder, &proc_grid);
  MPI_Comm_rank(proc_grid, &world_id);
  MPI_Cart_coords(proc_grid, world_id, ndim, coords);
  /*--------------------------------------------------------------------------*/

  /*Create a communicator for each row and column-----------------------------*/
  MPI_Comm_split(proc_grid, coords[0], coords[1], &proc_row);
  MPI_Comm_rank(proc_row, &row_id);
  MPI_Comm_split(proc_grid, coords[1], coords[0], &proc_col);
  MPI_Comm_rank(proc_col, &col_id);
  /*--------------------------------------------------------------------------*/

  /*Create a new datatype for the matrix partitioning-------------------------*/
  count = n_row/nodes_per_dim[0];
  blocklen = n_col/nodes_per_dim[1];
  stride = n_col;
  MPI_Type_vector(count, blocklen, stride, MPI_DOUBLE, &block);
  MPI_Type_commit(&block);
  /*--------------------------------------------------------------------------*/

  /*Create matrices and distribute block to sub matrices on each proc---------*/
  timer = MPI_Wtime();
  if (world_id == 0)
  {
    int rank;
    A = generate_matrix(n_row, n_col);
    B = generate_matrix(n_row, n_col);
    #ifdef REF
      Cref = reference_matrix(n_row, A, B);
    #endif

    for (i = 0; i < nodes_per_dim[0]; i++)
    {
      for (j = 0; j < nodes_per_dim[1]; j++)
      {
        pos[0] = i; pos[1] = j;
        idx = j*blocklen + (i*count*stride);
        MPI_Cart_rank(proc_grid, pos, &rank);
        MPI_Isend( &A[idx], 1, block, rank, 0, MPI_COMM_WORLD, &request);
        MPI_Isend( &B[idx], 1, block, rank, 1, MPI_COMM_WORLD, &request);
      }
    }
  }

  Asub = (double *)calloc(count*blocklen, sizeof(double));
  Bsub = (double *)calloc(count*blocklen, sizeof(double));
  Csub = (double *)calloc(count*blocklen, sizeof(double));
  MPI_Recv( Asub, count*blocklen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
  MPI_Recv( Bsub, count*blocklen, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
  MPI_Barrier(MPI_COMM_WORLD); //sync before deallocating A and B

  if (world_id == 0)
  {
    free(A);
    free(B);
  }
  /*--------------------------------------------------------------------------*/

  /* Initiate Fox's algorithm-------------------------------------------------*/
  Ainc = (double *)malloc(count*blocklen*sizeof(double));
  Binc = (double *)malloc(count*blocklen*sizeof(double));
  memcpy(Ainc, Asub, count*blocklen*sizeof(double));

  /* Figure out what target to send and receive from--------------------------*/
  int send_to = (col_id + nodes_per_dim[0] - 1)%nodes_per_dim[0];
  int recv_from = (col_id + 1)%nodes_per_dim[0];

  for (k = 0; k < nodes_per_dim[0]; k++)
  {
    /* PHASE1: Broadcast correct row block------------------------------------*/
    int m = 0;
    m = (col_id+k)%nodes_per_dim[0];
    MPI_Bcast(Ainc, count*blocklen, MPI_DOUBLE, m, proc_row);

    /* PHASE2 and 3: Shift Bsub and multiply----------------------------------*/
    if (k < nodes_per_dim[0]-1) {
      MPI_Isend( Bsub, count*blocklen, MPI_DOUBLE, send_to, 0, proc_col, &request);
      MPI_Recv( Binc, count*blocklen, MPI_DOUBLE, recv_from, 0, proc_col, &status);
    }

    mat_mul(count, Ainc, Bsub, Csub);
    MPI_Barrier(proc_col); // sync before shifting Bsub

    if (k < nodes_per_dim[0]-1) {
      memcpy(Ainc, Asub, count*blocklen*sizeof(double));
      memcpy(Bsub, Binc, count*blocklen*sizeof(double));
    }

  }
  /*--------------------------------------------------------------------------*/

  /* Send Csub on each proc to P0 and gather on P0----------------------------*/
  MPI_Isend( Csub, count*blocklen, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &request);

  if (world_id == 0)
  {
    C = (double *)calloc(n_row*n_col, sizeof(double));
    for (i = 0; i < nproc; i++) {
      MPI_Probe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
      MPI_Cart_coords(proc_grid, status.MPI_SOURCE, ndim, coords);
      MPI_Recv(&C[coords[1]*blocklen + (coords[0]*count*stride)], 1, block, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status );
    }
  }
  /*--------------------------------------------------------------------------*/
  MPI_Barrier(MPI_COMM_WORLD); //sync before deallocating C

  timer = MPI_Wtime() - timer;
  MPI_Reduce(&timer, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timer, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if (world_id == 0)
  {
    printf("length=%d\n", n_row);
    printf("processes=%d\n", nproc);
    printf("Max sort time=%f\n", max_time);
    printf("Avg. sort time=%f\n", avg_time/nproc);

    #ifdef REF
      int error_counter = check_matrix(n_col, Cref, C);
      printf("error_counter:%d\n", error_counter);
      free(Cref);
    #endif
    free(C);
  }

  free(Asub);
  free(Bsub);
  free(Csub);
  free(Ainc);
  free(Binc);
  MPI_Type_free(&block);
  MPI_Finalize();
  return 0;
}
