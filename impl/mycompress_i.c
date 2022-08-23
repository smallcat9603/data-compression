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

  int numtasks, rank, next, prev, tag=0;
  int numdata = 10;
  double* data = (double *)malloc(sizeof(double)*numdata);

  MPI_Request reqs[2];
  MPI_Status stats[2];

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  prev = rank - 1;
  next = rank + 1;
  if(rank == 0) prev = numtasks - 1;
  if(rank == numtasks - 1) next = 0;

  double send[numdata];
  for(int i=0; i<numdata; i++){
    if(rank == numtasks-1) send[i] = 0.0;
    else send[i] = rank + 1.0;
  }

  MPI_Irecv(data, numdata, MPI_DOUBLE, prev, tag, MPI_COMM_WORLD, &reqs[0]);
  MPI_Isend(send, numdata, MPI_DOUBLE, next, tag, MPI_COMM_WORLD, &reqs[1]);
  MPI_Waitall(2, reqs, stats);

  for(int i=0; i<numdata; i++){
    printf("Rank %d receives data[%d] from prev rank %d: %f\n", rank, i, prev, data[i]);
  }

  MPI_Finalize();

  free(data);

  return EXIT_SUCCESS;
}