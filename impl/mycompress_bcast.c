//
// my compress
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h> 
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "param.h"
#include "dataCompression.h"


int main(int argc, char** argv) {

  int num_procs, myrank;
  // double a, b;
  int num_data = 3;
  double* data = (double *)malloc(sizeof(double)*num_data);
  int tag = 0;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  // a = 0;
  // b = 0;

  if(myrank == 0) {
    // a = 1.0;
    double d = 1.0;
    for(int i = 0; i < num_data; i++)
    {
      data[i] = d;
      d += 1.0;
    }
  }
  // MPI_Bcast((void *)&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(data, num_data, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // printf("Process %d: a = %e, b = %e\n", myrank, a, b);
  for(int i = 0; i < num_data; i++)
  {
    printf("Process %d: data[%d] = %e\n", myrank, i, data[i]);
  }

  MPI_Finalize();

  free(data);

  return EXIT_SUCCESS;
}