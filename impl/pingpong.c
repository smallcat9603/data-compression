//
// Ping pong example with MPI_Send and MPI_Recv. Two processes ping pong a
// number back and forth, incrementing it until it reaches a given value.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.h"
#include "dataCompression.h"

struct vector
{
  float* p_data; //precise data
  char* c_data; //compressed data
};

int main(int argc, char** argv) {

  double start_time, end_time;

  const int PING_PONG_LIMIT = 10000;

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
  char output_filename[32];
  sprintf(output_filename, "%s%s", filename, suffix);
  FILE *fp = fopen(output_filename, "r");
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

  float* data_small = NULL;
  float min = toSmallDataset(data, &data_small, data_num);
  if (world_rank == 0)
  {
    printf("%f ", data[0]);
    printf("%f ", data_small[0]);
    printf("%f ", data[1]);
    printf("%f ", data_small[1]);
    printf("%f ", data[data_num-1]);
    printf("%f ", data_small[data_num-1]);
    printf("%f \n", min);
  }

  int array_float_len = myCompress(data, &array_float, &array_char, &array_char_displacement, data_num);

  struct vector msg; 
  int num_p = array_float_len, num_c = data_num-array_float_len;
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
      //printf("%d sent and incremented ping_pong_count %d to %d\n", world_rank, ping_pong_count, partner_rank);
      if(CT == 0)
      {
        MPI_Send(data, data_num, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD);
        //printf("%d sent data to %d\n", world_rank, partner_rank);
      }      
      if(CT == 1)
      {
        MPI_Send(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD);
        //printf("%d sent msg.p_data to %d\n", world_rank, partner_rank);
        MPI_Send(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD);
        //printf("%d sent msg.c_data to %d\n", world_rank, partner_rank);
      }   
    }
    else {
      MPI_Recv(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //printf("%d received ping_pong_count %d from %d\n", world_rank, ping_pong_count, partner_rank);
      if(CT == 0)
      {
        MPI_Recv(data, data_num, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("%d received data from %d\n", world_rank, partner_rank);
      }      
      if(CT == 1)
      {
        MPI_Recv(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("%d received msg.p_data from %d\n", world_rank, partner_rank);
        MPI_Recv(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("%d received msg.c_data from %d\n", world_rank, partner_rank);
      }
      if(ping_pong_count == PING_PONG_LIMIT && CT == 1)
      {
        float* decompressed_data = myDecompress(array_float, array_char, array_char_displacement, data_num);
        float gosa = 0;
        for(int i=0; i<data_num; i++)
        {
          gosa += fabs(decompressed_data[i]-data[i]);
        }
        gosa = gosa/data_num;
        printf("gosa = %f \n", gosa);
      }
    }
  }
  end_time = MPI_Wtime();

  if(world_rank == 0)
  {
    printf("rank = %d, elapsed = %f = %f - %f\n", world_rank, end_time-start_time, end_time, start_time);

    if (CT == 1)
    {
      float compress_ratio;
      compress_ratio = (float)(num_c*sizeof(char)+num_p*sizeof(float))/((num_c+num_p)*sizeof(float));
      printf("Compression rate (float, byte): %f \n", 1/compress_ratio); 
      compress_ratio = (float)(num_c*2+num_p*sizeof(float)*8)/((num_c+num_p)*sizeof(float)*8);
      printf("Compression rate (float, bit): %f \n", 1/compress_ratio); 
      compress_ratio = (float)(num_c*sizeof(char)+num_p*sizeof(double))/((num_c+num_p)*sizeof(double));
      printf("Compression rate (double, byte): %f \n", 1/compress_ratio); 
      compress_ratio = (float)(num_c*2+num_p*sizeof(double)*8)/((num_c+num_p)*sizeof(double)*8);
      printf("Compression rate (double, bit): %f \n", 1/compress_ratio);     

      compress_ratio = (3.0/(sizeof(float)*8))*((float)num_c/(num_c+num_p)) + calCompressRatio_bitwise_float(msg.p_data, num_p)*((float)num_p/(num_c+num_p));
      printf("Compression rate (bitwise, float): %f \n", 1/compress_ratio);        
      compress_ratio = (3.0/(sizeof(double)*8))*((float)num_c/(num_c+num_p)) + calCompressRatio_bitwise_double2(msg.p_data, num_p)*((float)num_p/(num_c+num_p));
      printf("Compression rate (bitwise, double): %f \n", 1/compress_ratio);    
    } 
  }

  MPI_Finalize();
}