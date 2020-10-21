#include <mpi.h>
#include <stdio.h>
#include <sys/time.h>

double timer() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int main(int argc, char *argv[]) {

  /*MPI_variables-------------------------------------------------------------*/
  int world_id, nproc, i, j;
  int coords[2], reorder = 1, ndim = 2, nodes_per_dim[2] = {0,0}, periods[2] = {0,0};
  MPI_Comm proc_grid;
  /*--------------------------------------------------------------------------*/

  /*Integral_variables--------------------------------------------------------*/
  const long int intervals = 100000L;
  int chunk_x, chunk_y;
  int istart, istop, jstart, jstop;
  double partial_sum, dx, dy, global_sum;
  double t = 0;
  /*--------------------------------------------------------------------------*/

  MPI_Init(&argc, &argv); /* Initialize MPI */

  /*Create a virtual 2D-grid topology-----------------------------------------*/
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Dims_create(nproc, ndim, nodes_per_dim);
  MPI_Cart_create(MPI_COMM_WORLD, ndim, nodes_per_dim, periods, reorder, &proc_grid);
  MPI_Comm_rank(proc_grid, &world_id);
  MPI_Cart_coords(proc_grid, world_id, ndim, coords);
  /*--------------------------------------------------------------------------*/

  /* Calculate how large a piece each process should do in x and y. */
  chunk_x  = intervals/nodes_per_dim[0];
  chunk_y  = intervals/nodes_per_dim[1];

  /* Calculate the start and stop index of each process in x and y. */
  istart = coords[0]*chunk_x+1;
  istop  = (coords[0] + 1)*chunk_x;
  jstart = coords[1]*chunk_y+1;
  jstop  = (coords[1] + 1)*chunk_y;
  //printf("(%d, %d) \t i[%d, %d] \t j[%d, %d]\n",coords[0], coords[1], istart, istop, jstart, jstop);

  /* Check that the entire integration range is covered */
  if (coords[0] == nodes_per_dim[0]-1 ) {
    istop = intervals;
  }
  if (coords[1] == nodes_per_dim[1]-1 ) {
    jstop = intervals;
  }

  /* Start timer if process has coordinates (0,0), ie the first process. */
  if (coords[0] == 0 && coords[1] == 0) {
    t = timer();
  }

  /* Compute the partial sum designated to this process */
  dx  = 1.0/intervals;
  dy  = 1.0/intervals;
  partial_sum = 0.0;
  for (i = istart; i <= istop; i++) {
    for(j = jstart; j <= jstop; ++j) {
      double x = dx*(i - 0.5);
      double y = dy*(j - 0.5);
      partial_sum += dx*dy*4.0/(1.0 + x*x + y*y);
    }
  }
  //printf("Partial sum = %.16f\n",  partial_sum);

  /* Reduce the partial sums to a global sum of all processes */
  global_sum = 0;
  MPI_Reduce(&partial_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  /* Stop the timer and print result */
  if (coords[0] == 0 && coords[1] == 0) {
    t = timer() - t;
    printf("Intervals=%ld \n",  intervals);
    printf("Integral=%.16f\n",  global_sum);
    printf("Time=%f s\n",  t);
    printf("Nodes_per_dim=(%d, %d) \n",  nodes_per_dim[0], nodes_per_dim[1]);
  }

  MPI_Finalize();

  return 0;
}
