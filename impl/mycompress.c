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
  double a, b;
  int tag = 0;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  a = 0;
  b = 0;

  if(myrank == 0)
  {
    a = 1.0;
    MPI_Send((void *)&a, 1, MPI_DOUBLE, 1, tag, MPI_COMM_WORLD);
  }
  else if(myrank == 1)
  {
    MPI_Recv((void *)&b, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
  }
  
  printf("Process %d: a = %e, b = %e\n", myrank, a, b);

  MPI_Finalize();
  
  return EXIT_SUCCESS;
}