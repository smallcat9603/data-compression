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

  double start_time, end_time;

  const int PING_PONG_LIMIT = 100;

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

  // read data file
  FILE *fp = fopen("testfloat_8_8_128.txt", "r");
  float *data = NULL; //data array
  int n; //data number = n-1
  for (n=0; !feof(fp); n++) 
  {
    data = (float *)(data?realloc(data,sizeof(float)*(n+1)):malloc(sizeof(float)));
    fscanf(fp, "%f", data+n);
    // printf("%f\t", data[n]);
  }
  fclose(fp);
  // free(data);

  // my compress
  float* array_float = NULL;
  char* array_char = NULL;
  int* array_char_displacement = NULL;
  int array_float_len = myCompress(data, &array_float, &array_char, &array_char_displacement);

  struct vector msg; 
  int num_p=array_float_len, num_c=data_num-array_float_len;
  // msg.p_data = (float*) malloc(sizeof(float)*num_p);
  // for (int i = 0; i < num_p; i++) {
  //   msg.p_data[i] = array_float[i];
  // }
  // msg.c_data = (char*) malloc(sizeof(char)*num_c);
  // for (int i = 0; i < num_c; i++) {
  //   msg.c_data[i] = array_char[i];
  // }
  msg.p_data = array_float;
  msg.c_data = array_char;

  int ping_pong_count = 0;
  int partner_rank = (world_rank + 1) % 2;
  start_time = MPI_Wtime();
  while (ping_pong_count < PING_PONG_LIMIT) {
    if (world_rank == ping_pong_count % 2) {
      // Increment the ping pong count before you send it
      ping_pong_count++;
      MPI_Send(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD);
      printf("%d sent and incremented ping_pong_count %d to %d\n", world_rank, ping_pong_count, partner_rank);
      // MPI_Send(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD);
      // printf("%d sent msg.p_data to %d\n", world_rank, partner_rank);
      // MPI_Send(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD);
      // printf("%d sent msg.c_data to %d\n", world_rank, partner_rank);
      MPI_Send(data, data_num, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD);
      printf("%d sent data to %d\n", world_rank, partner_rank);
    } else {
      MPI_Recv(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%d received ping_pong_count %d from %d\n", world_rank, ping_pong_count, partner_rank);
      // MPI_Recv(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // printf("%d received msg.p_data from %d\n", world_rank, partner_rank);
      // MPI_Recv(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // printf("%d received msg.c_data from %d\n", world_rank, partner_rank);
      MPI_Recv(data, data_num, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%d received data from %d\n", world_rank, partner_rank);
    }
  }
  end_time = MPI_Wtime();
  printf("rank = %d, elapsed = %f = %f - %f\n", world_rank, end_time-start_time, end_time, start_time);
  MPI_Finalize();
}