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
  double start_time_comp_byte, end_time_comp_byte, start_time_decomp_byte, end_time_decomp_byte;
  double start_time_comp_bit, end_time_comp_bit, start_time_decomp_bit, end_time_decomp_bit;
  double start_time_comp_bit_np, end_time_comp_bit_np, start_time_decomp_bit_np, end_time_decomp_bit_np;
  
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
  float min = toSmallDataset_float(data, &data_small, data_num);

  start_time_comp_byte = MPI_Wtime();
  int array_float_len = myCompress(data, &array_float, &array_char, &array_char_displacement, data_num);
  end_time_comp_byte = MPI_Wtime();

  float compress_ratio;

  float sz_comp_ratio = calcCompressionRatio_sz_float(data, data_num);
  float nolossy_performance = calcCompressionRatio_nolossy_performance_float(data, data_num);
  float nolossy_area = calcCompressionRatio_nolossy_area_float(data, data_num);
  printf("compression ratio: sz %f, nolossy_performance %f, nolossy_area %f \n", 1/sz_comp_ratio, 1/nolossy_performance, 1/nolossy_area);

  unsigned char* data_bits = NULL;
  //int flag = 0; //0, 1
  int bytes = 0; //total bytes of compressed data
  int pos = 8; //position of filled bit in last byte --> 87654321

  start_time_comp_bit = MPI_Wtime();
  myCompress_bitwise(data_small, data_num, &data_bits, &bytes, &pos);
  end_time_comp_bit = MPI_Wtime();
  //printf("test %d %d \n", bytes, pos);
  //printf("%.10f %.10f %.10f %.10f\n", data_small[0], data_small[1], data_small[2], data_small[data_num-1]);

  unsigned char* data_bits_np = NULL;
  //int flag = 0; //0, 1
  int bytes_np = 0; //total bytes of compressed data
  int pos_np = 8; //position of filled bit in last byte --> 87654321

  start_time_comp_bit_np = MPI_Wtime();
  myCompress_bitwise_np(data_small, data_num, &data_bits_np, &bytes_np, &pos_np);
  end_time_comp_bit_np = MPI_Wtime();  

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

  compress_ratio = (float)(num_c*sizeof(char)+num_p*sizeof(float))/((num_c+num_p)*sizeof(float));
  printf("Compression rate (float, byte): %f \n", 1/compress_ratio); 
  compress_ratio = (float)(num_c*2+num_p*sizeof(float)*8)/((num_c+num_p)*sizeof(float)*8);
  printf("Compression rate (float, bit): %f \n", 1/compress_ratio); 
  compress_ratio = (float)(num_c*sizeof(char)+num_p*sizeof(double))/((num_c+num_p)*sizeof(double));
  printf("Compression rate (double, byte): %f \n", 1/compress_ratio); 
  compress_ratio = (float)(num_c*2+num_p*sizeof(double)*8)/((num_c+num_p)*sizeof(double)*8);
  printf("Compression rate (double, bit): %f \n", 1/compress_ratio);     

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
      if(CT == 5)
      {
        MPI_Send(data_bits, bytes, MPI_CHAR, partner_rank, 4, MPI_COMM_WORLD);
      } 
      if(CT == 6)
      {
        MPI_Send(data_bits_np, bytes_np, MPI_CHAR, partner_rank, 5, MPI_COMM_WORLD);
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
      if(CT == 5)
      {
        MPI_Recv(data_bits, bytes, MPI_CHAR, partner_rank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      if(CT == 6)
      {
        MPI_Recv(data_bits_np, bytes_np, MPI_CHAR, partner_rank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }      
      if(ping_pong_count == PING_PONG_LIMIT)
      {
        if(CT == 1)
        {
          start_time_decomp_byte = MPI_Wtime();
          float* decompressed_data = myDecompress(array_float, array_char, array_char_displacement, data_num);
          end_time_decomp_byte = MPI_Wtime();

          float gosa = 0;
          for(int i=0; i<data_num; i++)
          {
            gosa += fabs(decompressed_data[i]-data[i]);
          }
          gosa = gosa/data_num;
          printf("gosa = %f \n", gosa);
        }
        if(CT == 5)
        {
          // for(int i=0; i<9; i++)
          // {
          //   for (int j=7; j>=0; j--)     //由低地址的位开始输出。
          //   {
          //     printf("%d", (data_bits[i] >> j) & 1);
          //   }
          //   printf(" ");          
          // }
          // printf("\n");
          start_time_decomp_bit = MPI_Wtime();
          float* decompressed_data = myDecompress_bitwise(data_bits, bytes, data_num);
          end_time_decomp_bit = MPI_Wtime();
          //printf("%.10f %.10f %.10f %.10f\n", decompressed_data[0], decompressed_data[1], decompressed_data[2], decompressed_data[data_num-1]);
          float gosa = 0;
          for(int i=0; i<data_num; i++)
          {
            gosa += fabs(decompressed_data[i]+min-data[i]);
          }
          gosa = gosa/data_num;
          printf("gosa = %f \n", gosa);          
        }
        if(CT == 6)
        {
          start_time_decomp_bit_np = MPI_Wtime();
          float* decompressed_data = myDecompress_bitwise_np(data_bits_np, bytes_np, data_num);
          end_time_decomp_bit_np = MPI_Wtime();
          //printf("%.10f %.10f %.10f %.10f\n", decompressed_data[0], decompressed_data[1], decompressed_data[2], decompressed_data[data_num-1]);
          float gosa = 0;
          for(int i=0; i<data_num; i++)
          {
            gosa += fabs(decompressed_data[i]+min-data[i]);
          }
          gosa = gosa/data_num;
          printf("gosa = %f \n", gosa);          
        }        
      }
    }
  }
  end_time = MPI_Wtime();

  if(world_rank == 0)
  {
    printf("rank = %d, elapsed = %f = %f - %f\n", world_rank, end_time-start_time, end_time, start_time);

    printf("Compression time (bytewise): %f \n", end_time_comp_byte-start_time_comp_byte); 
    printf("Compression time (bitwise): %f \n", end_time_comp_bit-start_time_comp_bit); 
    printf("Compression time (bitwise_np): %f \n", end_time_comp_bit_np-start_time_comp_bit_np); 

    if(CT == 1)
    {
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

      printf("Decompression time (bytewise): %f \n", end_time_decomp_byte-start_time_decomp_byte);  
    } 
    if(CT == 5)
    {
      compress_ratio = (float)(bytes*8)/(data_num*sizeof(float)*8);
      printf("Compression rate (bitwise, float): %f \n", 1/compress_ratio); 

      printf("Decompression time (bitwise): %f \n", end_time_decomp_bit-start_time_decomp_bit); 
    }
    if(CT == 6)
    {
      compress_ratio = (float)(bytes_np*8)/(data_num*sizeof(float)*8);
      printf("Compression rate (bitwise_np, float): %f \n", 1/compress_ratio); 

      printf("Decompression time (bitwise_np): %f \n", end_time_decomp_bit_np-start_time_decomp_bit_np); 
    }    
  }

  MPI_Finalize();
}