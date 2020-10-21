/**********************************************************************
 * Serial Quick sort
 * You can use this for Assignment 2
 *
 **********************************************************************/
#include <stdio.h>
#include <stdlib.h>

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


int main(int argc, char *argv[]) {
  double *data;
  int i, len;

  printf("Give number of elements: ");
  len = atoi(argv[1]);
  data = (double *)malloc(len*sizeof(double));

  // Generate random numbers
  for (i = 0; i < len; i++)
      data[i] = drand48();

  // Sort
  quicksort(data, 0, len-1);

  return 0;
}
