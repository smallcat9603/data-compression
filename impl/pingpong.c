//
// Ping pong example with MPI_Send and MPI_Recv. Two processes ping pong a
// number back and forth, incrementing it until it reaches a given value.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "param.h"
#include "dataCompression.h"

struct vector
{
  float* p_data; //precise data
  char* c_data; //compressed data
};

int main(int argc, char** argv) {

  const int PING_PONG_LIMIT = 10;

  // Initialize the MPI environment
  MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  if (world_size != 2) {
    fprintf(stderr, "World size must be two for %s\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  float* data = readFileFloat("testfloat_8_8_128.txt");

  struct vector msg; 
  int num_p=4, num_c=8;
  msg.p_data = (float*) malloc(sizeof(float)*num_p);
  msg.c_data = (char*) malloc(sizeof(char)*num_c);

  int ping_pong_count = 0;
  int partner_rank = (world_rank + 1) % 2;
  while (ping_pong_count < PING_PONG_LIMIT) {
    if (world_rank == ping_pong_count % 2) {
      // Increment the ping pong count before you send it
      ping_pong_count++;
      MPI_Send(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
      printf("%d sent and incremented ping_pong_count %d to %d\n", world_rank, ping_pong_count, partner_rank);
      for (int i = 0; i < num_p; i++) {
        msg.p_data[i] = i;
      }
      MPI_Send(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD);
      printf("%d sent msg.p_data to %d\n", world_rank, partner_rank);
      for (int i = 0; i < num_c; i++) {
        msg.c_data[i] = 'a';
      }
      MPI_Send(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD);
      printf("%d sent msg.c_data to %d\n", world_rank, partner_rank);
    } else {
      MPI_Recv(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%d received ping_pong_count %d from %d\n", world_rank, ping_pong_count, partner_rank);
      MPI_Recv(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%d received msg.p_data from %d\n", world_rank, partner_rank);
      MPI_Recv(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%d received msg.c_data from %d\n", world_rank, partner_rank);
    }
  }
  MPI_Finalize();
}