//
// Ping pong example with MPI_Send and MPI_Recv. Two processes ping pong a
// number back and forth, incrementing it until it reaches a given value.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "param.h"
#include "dataCompression.h"
// #include "/opt/sz191211/include/sz.h"
// #include "/opt/sz191211/include/rw.h"

struct vector
{
  float* p_data; //precise data
  //double* p_data;  //switch to double
  char* c_data; //compressed data
};

int main(int argc, char** argv) {

  double start_time, end_time;
  double start_time_comp_byte, end_time_comp_byte, start_time_decomp_byte, end_time_decomp_byte;
  double start_time_comp_bit, end_time_comp_bit, start_time_decomp_bit, end_time_decomp_bit;
  double start_time_comp_bit_np, end_time_comp_bit_np, start_time_decomp_bit_np, end_time_decomp_bit_np;
  double start_time_comp_sz, end_time_comp_sz, start_time_decomp_sz, end_time_decomp_sz;
  double start_time_comp_bit_mask, end_time_comp_bit_mask, start_time_decomp_bit_mask, end_time_decomp_bit_mask;
  
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
  //double *data = NULL; //switch to double
  int n; //data number = n-1
  for (n=0; !feof(fp); n++) 
  {
    data = (float *)(data?realloc(data,sizeof(float)*(n+1)):malloc(sizeof(float)));
    //data = (double *)(data?realloc(data,sizeof(double)*(n+1)):malloc(sizeof(double))); //switch to double
    fscanf(fp, "%f", data+n);
    // printf("%f\t", data[n]);
  }
  fclose(fp);
  // free(data);
  int data_num = n - 1;

  float compress_ratio;

  // float sz_comp_ratio = calcCompressionRatio_sz_float(data, data_num);
  // float nolossy_performance = calcCompressionRatio_nolossy_performance_float(data, data_num);
  // float nolossy_area = calcCompressionRatio_nolossy_area_float(data, data_num);
  // printf("compression ratio: sz %f, nolossy_performance %f, nolossy_area %f \n", 1/sz_comp_ratio, 1/nolossy_performance, 1/nolossy_area);  

  //sz compression
  start_time_comp_sz = MPI_Wtime();
  char* binfile = filename bin_suffix; 
  writetobinary_float(binfile, data, data_num); //.txt --> .dat
  //writetobinary_double(binfile, data, data_num); //switch to double
  char sz_comp_cmd[64];
  sprintf(sz_comp_cmd, "%s%g%s%s%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, filename, sz_comp_cmd_suffix2, data_num);
  //sprintf(sz_comp_cmd, "%s%g%s%s%s%d", sz_comp_cmd_prefix_double, absErrorBound, sz_comp_cmd_suffix1, filename, sz_comp_cmd_suffix2, data_num); //switch to double
  //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
  int iret_comp = system(sz_comp_cmd); //.dat --> .dat.sz
  char* binfile_sz = filename bin_suffix sz_suffix;
  int bytes_sz = 0;
  unsigned char* data_bits_sz = readfrombinary_char(binfile_sz, &bytes_sz);
  end_time_comp_sz = MPI_Wtime();

  // my compress bytewise
  float* array_float = NULL;
  //double* array_double = NULL; //switch to double
  char* array_char = NULL;
  int* array_char_displacement = NULL;

  start_time_comp_byte = MPI_Wtime();
  int array_float_len = myCompress(data, &array_float, &array_char, &array_char_displacement, data_num);
  //int array_double_len = myCompress_double(data, &array_double, &array_char, &array_char_displacement, data_num); //switch to double
  end_time_comp_byte = MPI_Wtime();

  struct vector msg; 
  int num_p = array_float_len, num_c = data_num-array_float_len;
  //int num_p = array_double_len, num_c = data_num-array_double_len; //switch to double

  // msg.p_data = (float*) malloc(sizeof(float)*num_p);
  // for (int i = 0; i < num_p; i++) {
  //   msg.p_data[i] = array_float[i];
  // }
  // msg.c_data = (char*) malloc(sizeof(char)*num_c);
  // for (int i = 0; i < num_c; i++) {
  //   msg.c_data[i] = array_char[i];
  // }

  msg.p_data = array_float;
  //msg.p_data = array_double; //switch to double
  msg.c_data = array_char;

  // compress_ratio = (float)(num_c*sizeof(char)+num_p*sizeof(float))/((num_c+num_p)*sizeof(float));
  // printf("Compression rate (float, byte): %f \n", 1/compress_ratio); 
  // compress_ratio = (float)(num_c*2+num_p*sizeof(float)*8)/((num_c+num_p)*sizeof(float)*8);
  // printf("Compression rate (float, bit): %f \n", 1/compress_ratio); 
  // compress_ratio = (float)(num_c*sizeof(char)+num_p*sizeof(double))/((num_c+num_p)*sizeof(double));
  // printf("Compression rate (double, byte): %f \n", 1/compress_ratio); 
  // compress_ratio = (float)(num_c*2+num_p*sizeof(double)*8)/((num_c+num_p)*sizeof(double)*8);
  // printf("Compression rate (double, bit): %f \n", 1/compress_ratio);    

  // my compress bitwise
  float* data_small = NULL;
  //double* data_small = NULL; //switch to double
  float min = toSmallDataset_float(data, &data_small, data_num);
  //double min = toSmallDataset_double(data, &data_small, data_num); //switch to double

  // fp = fopen("bit.txt", "w");
  // for(int i=0; i<data_num; i++)
  // {
  //   float float10 = data[i] - min;
  //   char float_arr[32+1];
  //   floattostr(&float10, float_arr);
  //   for(int j=0; j<32; j++)
  //   {
  //     fprintf(fp, "%d ", float_arr[j]-'0');
  //   }
  //   fprintf(fp, "\n");
  // }
  // fclose(fp);  

  unsigned char* data_bits = NULL;
  //int flag = 0; //0, 1
  int bytes = 0; //total bytes of compressed data
  int pos = 8; //position of filled bit in last byte --> 87654321

  start_time_comp_bit = MPI_Wtime();
  myCompress_bitwise(data_small, data_num, &data_bits, &bytes, &pos);
  //myCompress_bitwise_double(data_small, data_num, &data_bits, &bytes, &pos); //switch to double
  end_time_comp_bit = MPI_Wtime();
  //printf("test %d %d \n", bytes, pos);
  //printf("%.10f %.10f %.10f %.10f\n", data_small[0], data_small[1], data_small[2], data_small[data_num-1]);

  // my compress bitwise with no prediction
  unsigned char* data_bits_np = NULL;
  //int flag = 0; //0, 1
  int bytes_np = 0; //total bytes of compressed data
  int pos_np = 8; //position of filled bit in last byte --> 87654321

  start_time_comp_bit_np = MPI_Wtime();
  myCompress_bitwise_np(data_small, data_num, &data_bits_np, &bytes_np, &pos_np);
  //myCompress_bitwise_double_np(data_small, data_num, &data_bits_np, &bytes_np, &pos_np); //switch to double
  end_time_comp_bit_np = MPI_Wtime();    

  // my compress bitmask-based bitwise
  unsigned char* data_bits_mask = NULL;
  //int flag = 0; //0, 1
  int bytes_mask = 0; //total bytes of compressed data
  int pos_mask = 8; //position of filled bit in last byte --> 87654321
  int type = 0;

  float medium = med_dataset_float(data_small, data_num, &type);
  //double medium = med_dataset_double(data_small, data_num, &type); //switch to double
  // printf("medium = %f\n", medium);
  // printf("type = %d\n", type);
  char float_arr[32+1];
  //char double_arr[64+1]; //switch to double
  floattostr(&medium, float_arr);
  //doubletostr(&medium, double_arr); //switch to double
  char mask[1+8+8];
  //char mask[1+11+8]; //switch to double
  strncpy(mask, float_arr, 1+8+8);
  //strncpy(mask, double_arr, 1+11+8); //switch to double
  // for(int n=0; n<17; n++)
  // {
  //   printf("%c", mask[n]);
  // }
  // printf("\n");
  start_time_comp_bit_mask = MPI_Wtime();
  myCompress_bitwise_mask(data_small, data_num, &data_bits_mask, &bytes_mask, &pos_mask, type, mask);
  //myCompress_bitwise_double_mask(data_small, data_num, &data_bits_mask, &bytes_mask, &pos_mask, type, mask); //switch to double
  // printf("bytes_mask = %d, improvement = %f\n", bytes_mask, (float)bytes_mask/bytes);
  end_time_comp_bit_mask = MPI_Wtime();     

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
        //MPI_Send(data, data_num, MPI_DOUBLE, partner_rank, 3, MPI_COMM_WORLD); //switch to double
        //printf("%d sent data to %d\n", world_rank, partner_rank);
      }      
      if(CT == 1)
      {
        MPI_Send(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD);
        //MPI_Send(msg.p_data, num_p, MPI_DOUBLE, partner_rank, 1, MPI_COMM_WORLD); //switch to double
        //printf("%d sent msg.p_data to %d\n", world_rank, partner_rank);
        MPI_Send(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD);
        //printf("%d sent msg.c_data to %d\n", world_rank, partner_rank);
      }
      if(CT == 4)
      {
        MPI_Send(data_bits_sz, bytes_sz, MPI_UNSIGNED_CHAR, partner_rank, 6, MPI_COMM_WORLD);
      }       
      if(CT == 5)
      {
        MPI_Send(data_bits, bytes, MPI_CHAR, partner_rank, 4, MPI_COMM_WORLD);
      } 
      if(CT == 6)
      {
        MPI_Send(data_bits_np, bytes_np, MPI_CHAR, partner_rank, 5, MPI_COMM_WORLD);
      }        
      if(CT == 7)
      {
        MPI_Send(data_bits_mask, bytes_mask, MPI_CHAR, partner_rank, 7, MPI_COMM_WORLD);
      }   
    }
    else {
      MPI_Recv(&ping_pong_count, 1, MPI_INT, partner_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //printf("%d received ping_pong_count %d from %d\n", world_rank, ping_pong_count, partner_rank);
      if(CT == 0)
      {
        MPI_Recv(data, data_num, MPI_FLOAT, partner_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Recv(data, data_num, MPI_DOUBLE, partner_rank, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //switch to double
        //printf("%d received data from %d\n", world_rank, partner_rank);
      }      
      if(CT == 1)
      {
        MPI_Recv(msg.p_data, num_p, MPI_FLOAT, partner_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //MPI_Recv(msg.p_data, num_p, MPI_DOUBLE, partner_rank, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); //switch to double
        //printf("%d received msg.p_data from %d\n", world_rank, partner_rank);
        MPI_Recv(msg.c_data, num_c, MPI_CHAR, partner_rank, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("%d received msg.c_data from %d\n", world_rank, partner_rank);
      }
      if(CT == 4)
      {
        MPI_Recv(data_bits_sz, bytes_sz, MPI_UNSIGNED_CHAR, partner_rank, 6, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }        
      if(CT == 5)
      {
        MPI_Recv(data_bits, bytes, MPI_CHAR, partner_rank, 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
      if(CT == 6)
      {
        MPI_Recv(data_bits_np, bytes_np, MPI_CHAR, partner_rank, 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }    
      if(CT == 7)
      {
        MPI_Recv(data_bits_mask, bytes_mask, MPI_CHAR, partner_rank, 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }           
      if(ping_pong_count == PING_PONG_LIMIT)
      {
        if(CT == 1)
        {
          start_time_decomp_byte = MPI_Wtime();
          float* decompressed_data = myDecompress(array_float, array_char, array_char_displacement, data_num);
          //double* decompressed_data = myDecompress_double(array_double, array_char, array_char_displacement, data_num); //switch to double
          end_time_decomp_byte = MPI_Wtime();

          float gosa = 0;
          for(int i=0; i<data_num; i++)
          {
            gosa += fabs(decompressed_data[i]-data[i]);
          }
          gosa = gosa/data_num;
          printf("gosa = %f \n", gosa);
        }
        if(CT == 4)
        {
          start_time_decomp_sz = MPI_Wtime();
          char* binfile_zs = filename bin_suffix zs_suffix;
          writetobinary_char(binfile_zs, data_bits_sz, bytes_sz); //.dat.zs
          char sz_decomp_cmd[64];
          sprintf(sz_decomp_cmd, "%s%s%s%d", sz_decomp_cmd_prefix, filename, sz_decomp_cmd_suffix, data_num);
          //sprintf(sz_decomp_cmd, "%s%s%s%d", sz_decomp_cmd_prefix_double, filename, sz_decomp_cmd_suffix, data_num); //switch to double
          //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
          int iret_decomp = system(sz_decomp_cmd); //.dat.zs --> .dat.zs.out
          char* binfile_out = filename bin_suffix zs_suffix out_suffix;
          char* txtfile = filename bin_suffix zs_suffix out_suffix suffix;  
          float* decompressed_data = readfrombinary_writetotxt_float(binfile_out, txtfile, data_num);
          //double* decompressed_data = readfrombinary_writetotxt_double(binfile_out, txtfile, data_num); //switch to double
          end_time_decomp_sz = MPI_Wtime();

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
          start_time_decomp_bit = MPI_Wtime();
          float* decompressed_data = myDecompress_bitwise(data_bits, bytes, data_num);
          //double* decompressed_data = myDecompress_bitwise_double(data_bits, bytes, data_num); //switch to double
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
          //double* decompressed_data = myDecompress_bitwise_double_np(data_bits_np, bytes_np, data_num); //switch to double
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
        if(CT == 7)
        {
          start_time_decomp_bit_mask = MPI_Wtime();
          float* decompressed_data = myDecompress_bitwise_mask(data_bits_mask, bytes_mask, data_num, type, mask);
          //double* decompressed_data = myDecompress_bitwise_double_mask(data_bits_mask, bytes_mask, data_num, type, mask); //switch to double
          end_time_decomp_bit_mask = MPI_Wtime();
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
    printf("Compression time (sz): %f \n", end_time_comp_sz-start_time_comp_sz); 
    printf("Compression time (bitwise_mask): %f \n", end_time_comp_bit_mask-start_time_comp_bit_mask); 

    if(CT == 1)
    {    
      printf("Decompression time (bytewise): %f \n", end_time_decomp_byte-start_time_decomp_byte);  
      // compress_ratio = (3.0/(sizeof(float)*8))*((float)num_c/(num_c+num_p)) + calCompressRatio_bitwise_float(msg.p_data, num_p)*((float)num_p/(num_c+num_p));
      // printf("Compression rate (bitwise, float): %f \n", 1/compress_ratio);        
      // compress_ratio = (3.0/(sizeof(double)*8))*((float)num_c/(num_c+num_p)) + calCompressRatio_bitwise_double2(msg.p_data, num_p)*((float)num_p/(num_c+num_p));
      // printf("Compression rate (bitwise, double): %f \n", 1/compress_ratio); 
    } 
    if(CT == 4)
    {
      printf("Decompression time (sz): %f \n", end_time_decomp_sz-start_time_decomp_sz); 
      compress_ratio = (float)(bytes_sz*8)/(data_num*sizeof(float)*8);
      //compress_ratio = (float)(bytes_sz*8)/(data_num*sizeof(double)*8); //switch to double
      printf("Compression rate (sz): %f \n", 1/compress_ratio); 
    }
    if(CT == 5)
    {
      printf("Decompression time (bitwise): %f \n", end_time_decomp_bit-start_time_decomp_bit); 
      compress_ratio = (float)(bytes*8)/(data_num*sizeof(float)*8);
      //compress_ratio = (float)(bytes*8)/(data_num*sizeof(double)*8); //switch to double
      printf("Compression rate (bitwise): %f \n", 1/compress_ratio); 
    }
    if(CT == 6)
    {
      printf("Decompression time (bitwise_np): %f \n", end_time_decomp_bit_np-start_time_decomp_bit_np); 
      compress_ratio = (float)(bytes_np*8)/(data_num*sizeof(float)*8);
      //compress_ratio = (float)(bytes_np*8)/(data_num*sizeof(double)*8); //switch to double
      printf("Compression rate (bitwise_np): %f \n", 1/compress_ratio); 
    }    
    if(CT == 7)
    {
      printf("Decompression time (bitwise_mask): %f \n", end_time_decomp_bit_mask-start_time_decomp_bit_mask); 
      compress_ratio = (float)(bytes_mask*8)/(data_num*sizeof(float)*8);
      //compress_ratio = (float)(bytes_mask*8)/(data_num*sizeof(double)*8); //switch to double
      printf("Compression rate (bitwise_mask): %f (improvement = %f) \n", 1/compress_ratio, (float)bytes/bytes_mask); 
    }     
  }

  MPI_Finalize();
}