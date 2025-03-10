/**********************************************************************
 * Parallel Quick sort
 * Based on serial quicksort provided by Jarmo Rantakokko, which includes
 * the functions swap() and quicksort().
 * The functions merge_array and find_pivot are based on boilerplate code
 * for the algorithms mergesort and binary search, both found online.
 **********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

//#define REFERENCE

void merge_arrays(double* v1, double* v2, int start1, int stop1, int start2, int stop2, double* merge)
{
  int i = start1, j = start2, k = 0;

  while (i < stop1 && j < stop2)
  {
      if ( *(v1+i) < *(v2+j) ) {
        *(merge+k) = *(v1+i);
        k++; i++;
      }
      else
      {
        *(merge+k) = *(v2+j);
        k++; j++;
      }

  }

  while (i < stop1)
  {
    *(merge+k) = *(v1+i);
    k++; i++;
  }

  while (j < stop2)
  {
    *(merge+k) = *(v2+j);
    k++; j++;
  }

}


int find_pivot(double* v, int len, double value)
{
  int mid = (len-1)/2;
  int first = 0;
  int last = len-1;

  while (first <= last)
  {
    if (*(v+mid) <= value && *(v+mid+1) >= value)
    {
      return (mid + 1);
    }
    else
    {
      if ( *(v+mid) < value)
      {
        first = mid + 1;
      }
      else if ( *(v+mid) == value)
      {
         return mid;
      }
      else
      {
        last = mid - 1;
      }

      mid = (first + last)/2;
    }
  }
   if (first > last) {
     return len;
   }

   return 0;
 }


void swap(double *v, int i, int j)
{
  double temp;
  temp = v[i];
  v[i] = v[j];
  v[j] = temp;
}


void quicksort(double *v, int left, int right)
{
  int i, prev_index;
  if(left >= right)
    return;
  swap(v, left, (left+right)/2);
  prev_index = left;

  for(i = left+1; i <= right; i++)
    if(v[i] < v[left])
    {
      swap(v, ++prev_index, i);
    }

  swap(v, left, prev_index);

  quicksort(v, left, prev_index-1);
  quicksort(v, prev_index+1, right);
}


int counter = 0;
double* p_quicksort( double* data, int elements, MPI_Comm comm)
{
  counter++;
  int i, pid, nproc, pivot_elem, swap_target, inc_elems, merge_size;
  MPI_Comm_rank(comm, &pid);
  MPI_Comm_size(comm, &nproc);
  double local_median, *merge, *in, *final_data_buffer;
  double medians[nproc];
  MPI_Status status;
  MPI_Request request;
  MPI_Comm split_comm;

  /* CHECK IF CURRENT ITERATION ONLY CONTAINS ONE PROCESS */
  if (nproc == 1)
  {
    MPI_Isend( data, elements, MPI_DOUBLE, pid, 0, MPI_COMM_WORLD, &request);
    return data;
  }
  else
  {
    /* GENERATE PIVOT, MEDIAN OF MEDIANS */
    local_median = data[elements/2];
    MPI_Gather(&local_median, 1, MPI_DOUBLE, &medians[0], 1, MPI_DOUBLE, 0, comm);

    if (pid == 0)
    {
      quicksort(medians, 0, nproc-1);
      local_median = medians[nproc/2];
    }
    MPI_Bcast(&local_median, 1, MPI_DOUBLE, 0, comm);
    pivot_elem = find_pivot(data, elements, local_median);

    /* GENERATE PIVOT, MEDIAN OF P0 */
    // if (pid == 0) {
    //   local_median = data[elements/2];
    // }
    // MPI_Bcast(&local_median, 1, MPI_DOUBLE, 0, comm);
    // pivot_elem = find_pivot(data, elements, local_median);

    /* --------------------------------LOW----------------------------------- */
    if (pid < nproc/2)
    {
      /* RECEIVE AND SEND */
      swap_target = pid+nproc/2;
      MPI_Probe(swap_target, 0, comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &inc_elems);

      in = (double *)malloc( inc_elems*sizeof(double) );
      MPI_Recv( in, inc_elems, MPI_DOUBLE, swap_target, 0, comm, &status);

      MPI_Isend( (data+pivot_elem), elements-pivot_elem, MPI_DOUBLE, swap_target, 0, comm, &request);
      MPI_Wait(&request, &status);
      //printf("COUNT[%d]:PID[%d]:LOW:SOURCE=%d \t inc_elems=%d \t local_median=%f \t pivot_elem=%d \n", counter, pid, status.MPI_SOURCE, inc_elems, local_median, pivot_elem);

      /* MERGE */
      merge_size = pivot_elem + inc_elems;
      merge = (double *)malloc( merge_size*sizeof(double) );
      merge_arrays(in, data, 0, inc_elems, 0, pivot_elem, merge);
      free(data);
      free(in);

      /* SPLIT COMM */
      MPI_Comm_split(comm, 0, pid, &split_comm);

      /* RECURSION */
      final_data_buffer = p_quicksort( merge, merge_size, split_comm);
      MPI_Comm_free(&split_comm);
    }
    /* --------------------------------HIGH--------------------------------- */
    if (pid >= nproc/2)
    {
      /* SEND AND RECEIVE */
      swap_target = pid-nproc/2;
      MPI_Isend( data, pivot_elem, MPI_DOUBLE, swap_target, 0, comm, &request);
      MPI_Wait(&request, &status);

      MPI_Probe(swap_target, 0, comm, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &inc_elems);

      in = (double *)malloc( inc_elems*sizeof(double) );
      MPI_Recv( in, inc_elems, MPI_DOUBLE, swap_target, 0, comm, &status);
      //printf("COUNT[%d]:PID[%d]:HIGH:SOURCE=%d\t elements=%d \t inc_elems=%d \t local_median=%f \t pivot_elem=%d \n", counter, pid, status.MPI_SOURCE, elements, inc_elems, local_median, pivot_elem);

      /* MERGE */
      merge_size = elements - pivot_elem + inc_elems;
      merge = (double *)malloc( merge_size*sizeof(double) );
      merge_arrays(in, data, 0, inc_elems, pivot_elem, elements, merge);
      free(data);
      free(in);

      /* SPLIT COMM */
      MPI_Comm_split(comm, 1, pid, &split_comm);

      /* RECURSION */
      final_data_buffer = p_quicksort( merge, merge_size, split_comm);
      MPI_Comm_free(&split_comm);
    }
  }
  return final_data_buffer;
}


/*----------------------------------------------------------------------------*/
/*********************************MAIN*****************************************/
/*----------------------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  if (argc != 2)
  {
    printf("Invalid parameter input!\n");
    printf("Input parameters on the form:\n");
    printf("./p_qsort array_length \n");
    return 1;
    exit(1);
  }
  double *global, *data, *reference, timer, avg_time, max_time;
  int i, pid, nproc, elements, inc_elems, remain, *element_array, *offset_array;
  int len = atoi(argv[1]);
  MPI_Status status;
  MPI_Request request;

  MPI_Init(&argc, &argv);
  timer = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  elements = len/nproc;
  remain = len%nproc;

  /* Setup a local data buffer for each process */
  data = (double *)malloc( elements * sizeof(double) );

  /* Setup a global array for the final sort */
  if (pid == 0)
  {
    global = (double *)malloc(len*sizeof(double));
    for (i = 0; i < len; i++)
    {
      global[i] = drand48();
    }
  }

  /* HANDLE INPUT WHERE len%nproc != 0 */
  /* Get the remainder, distribute it among the processors and send the
  number of elements for each process to  P0 */
  if (remain > pid)
  {
    elements = elements+1;
  }
  MPI_Isend( &elements, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &request);

  /* P0 receives info about how many elements each process should receive in
  the data[] buffer. This and the offset in global[] is stored in a separate
  array to be used with MPI_Scatterv(). */
  if (pid == 0)
  {
    int offset = 0;
    element_array = (int *)malloc(nproc*sizeof(int));
    offset_array = (int *)malloc(nproc*sizeof(int));
    for ( i = 0; i < nproc; i++)
    {
      MPI_Irecv( element_array+i, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);
      *(offset_array+i) = offset;
      offset += *(element_array+i);
    }
  }

  /* Distribute the global array to all processes */
  MPI_Scatterv(global, element_array, offset_array, MPI_DOUBLE, data, elements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (pid == 0)
  {
    free(global);
    free(element_array);
    free(offset_array);
  }
  quicksort(data, 0, elements-1);
  data = p_quicksort(data, elements, MPI_COMM_WORLD);

  /* receive sorted local data[] into global[] on P0 */
  if (pid == 0)
  {
    global = (double *)malloc(len*sizeof(double));
    int offset = 0;
    for ( i = 0; i < nproc; i++)
    {
      MPI_Probe(i, 0, MPI_COMM_WORLD, &status);
      MPI_Get_count(&status, MPI_DOUBLE, &inc_elems);
      MPI_Irecv( global+offset, inc_elems, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);
      offset += inc_elems;
    }
  }

  timer = MPI_Wtime() - timer;
  MPI_Reduce(&timer, &avg_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timer, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  free(data);

  if (pid == 0)
  {
    printf("PID[%d]:Max sort time=%f\n", pid, max_time);
    printf("PID[%d]:Avg. sort time=%f\n", pid, avg_time/nproc);
  }

  /* Generate reference[] and compare received data in globa[] with reference[] */
  if ( pid==0 )
  {
    #ifdef REFERENCE
    int counter = 0;
    reference = (double *)malloc(len*sizeof(double));
    for (i = 0; i < len; i++)
    {
      reference[i] = global[i];
    }
    quicksort(reference, 0, len-1);

    for (i = 0; i < len; i++)
    {
      if (global[i] != reference[i])
      {
        counter++;
      }
    }
    printf("Faulty elements=%d\n", counter);
    free(reference);
    #endif

    free(global);
  }

  MPI_Finalize();
  return 0;
}
