#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  char input_txt[64], input_bi[128], output_bc[128], output_txt[128]; 

  if(argc < 2)
  {
    printf("Test case: ./mycompress_bitop_float [input_txt]\n");
    printf("Example: ./mycompress_bitop_float input.txt\n");
    exit(0);
  }
  
  sprintf(input_txt, "%s", argv[1]);
  printf("input_txt = %s\n", input_txt); 
  sprintf(input_bi, "%s.bi", input_txt);
  sprintf(output_bc, "%s.bc", input_txt);
  sprintf(output_txt, "%s.bop.txt", input_txt);

  FILE *fp = fopen(input_txt, "r");
  float *data = NULL; //data array
  int n; //data number = n-1
  for (n=0; !feof(fp); n++) 
  {
    data = (float *)(data?realloc(data,sizeof(float)*(n+1)):malloc(sizeof(float)));
    fscanf(fp, "%f", data+n);
  }
  fclose(fp);
  int num = n-1;

  writetobinary_float(input_bi, data, num); //.txt --> .bi

  float* data_small = NULL;
  float min = toSmallDataset_float(data, &data_small, num);

  unsigned char* data_bits_op = NULL;
  int bytes_op = 0; //total bytes of compressed data
  int pos_op = 8; //position of filled bit in last byte --> 87654321

  double start_time_comp = clock();
  myCompress_bitwise_op(data_small, num, &data_bits_op, &bytes_op, &pos_op);
  double end_time_comp = clock();

  writetobinary_char(output_bc, data_bits_op, bytes_op); //.bc

  int bytes_decomp = 0;
  unsigned char* data_bits_decomp = readfrombinary_char(output_bc, &bytes_decomp);

  double start_time_decomp = clock();
  float* decompressed_data = myDecompress_bitwise_op(data_bits_decomp, bytes_decomp, num);
  double end_time_decomp = clock();

  fp = fopen(output_txt, "w");
  for(int i=0; i<num; i++)
  {
    fprintf(fp, "%f\n", decompressed_data[i]+min);
  }  
  fclose(fp);
  printf("%sに保存しました。\n", output_txt);

  printf("absErrorBound: %f \n", absErrorBound); 

  float compress_ratio = (float)(bytes_op*8)/(num*sizeof(float)*8);
  printf("Compression rate: %f \n", 1/compress_ratio); 
  printf("Compression time: %f sec \n", (end_time_comp-start_time_comp)/CLOCKS_PER_SEC); 
  printf("Decompression time: %f sec \n", (end_time_decomp-start_time_decomp)/CLOCKS_PER_SEC); 

  printf("done\n");

  free(data);
  
  return 0;

}