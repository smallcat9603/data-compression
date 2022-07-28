#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  char input_txt[64], input_bi[128], output_bc[128], output_txt[128]; 

  if(argc < 2)
  {
    printf("Test case: ./mycompress_bitmask_double [input_txt]\n");
    printf("Example: ./mycompress_bitmask_double input.txt\n");
    exit(0);
  }
  
  sprintf(input_txt, "%s", argv[1]);
  printf("input_txt = %s\n", input_txt); 
  sprintf(input_bi, "%s.bi", input_txt);
  sprintf(output_bc, "%s.bc", input_txt);
  sprintf(output_txt, "%s.bm.txt", input_txt);

  FILE *fp = fopen(input_txt, "r");
  double *data = NULL; //data array
  int n; //data number = n-1
  for (n=0; !feof(fp); n++) 
  {
    data = (double *)(data?realloc(data,sizeof(double)*(n+1)):malloc(sizeof(double)));
    fscanf(fp, "%lf", data+n);
  }
  fclose(fp);
  int num = n-1;

  writetobinary_double(input_bi, data, num); //.txt --> .bi

  double* data_small = NULL;
  double min = toSmallDataset_double(data, &data_small, num);

  unsigned char* data_bits_mask = NULL;
  int bytes_mask = 0; //total bytes of compressed data
  int pos_mask = 8; //position of filled bit in last byte --> 87654321
  int type = 0;

  double medium = med_dataset_double(data_small, num, &type);
  char double_arr[64+1];
  doubletostr(&medium, double_arr);
  char mask[1+11+8];
  strncpy(mask, double_arr, 1+11+8); 

  double start_time_comp = clock();
  myCompress_bitwise_double_mask(data_small, num, &data_bits_mask, &bytes_mask, &pos_mask, type, mask);
  double end_time_comp = clock();

  writetobinary_char(output_bc, data_bits_mask, bytes_mask); //.bc

  int bytes_mask_decomp = 0;
  unsigned char* data_bits_mask_decomp = readfrombinary_char(output_bc, &bytes_mask_decomp);

  double start_time_decomp = clock();
  double* decompressed_data = myDecompress_bitwise_double_mask(data_bits_mask_decomp, bytes_mask_decomp, num, type, mask);
  double end_time_decomp = clock();

  fp = fopen(output_txt, "w");
  for(int i=0; i<num; i++)
  {
    fprintf(fp, "%f\n", decompressed_data[i]+min);
  }  
  fclose(fp);
  printf("%sに保存しました。\n", output_txt);

  printf("absErrorBound: %f \n", absErrorBound); 

  float compress_ratio = (float)(bytes_mask*8)/(num*sizeof(double)*8);
  printf("Compression rate: %f \n", 1/compress_ratio); 
  printf("Compression time: %f sec \n", (end_time_comp-start_time_comp)/CLOCKS_PER_SEC); 
  printf("Decompression time: %f sec \n", (end_time_decomp-start_time_decomp)/CLOCKS_PER_SEC); 

  printf("done\n");

  free(data);
  
  return 0;

}