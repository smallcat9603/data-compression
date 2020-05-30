#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  char input[64], output[64]; 
  char binfile[64], decomp[64];
  
  if(argc < 1)
  {
  printf("Test case: mycompress [srcFilePath]\n");
  printf("Example: mycompress test.txt\n");
  exit(0);
  }
  
  sprintf(input, "%s", argv[1]);
  
  printf("input=%s\n", input); 
  sprintf(output, "%s.bc", input);
  sprintf(binfile, "%s.bi", input);

  FILE *fp = fopen(input, "r");
  float *data = NULL; //data array
  int n; //data number = n-1
  for (n=0; !feof(fp); n++) 
  {
    data = (float *)(data?realloc(data,sizeof(float)*(n+1)):malloc(sizeof(float)));
    fscanf(fp, "%f", data+n);
  }
  fclose(fp);
  int data_num = n-1;

  writetobinary_float(binfile, data, data_num); //.txt --> .dat

  float* data_small = NULL;
  float min = toSmallDataset_float(data, &data_small, data_num);

  float compress_ratio;

  unsigned char* data_bits_comp = NULL;
  int bytes_comp = 0; //total bytes of compressed data
  int pos = 8; //position of filled bit in last byte --> 87654321
  myCompress_bitwise(data_small, data_num, &data_bits_comp, &bytes_comp, &pos);

  int bytes_decomp = 0;
  unsigned char* data_bits_decomp = readfrombinary_char(output, &bytes_decomp);
  
  writetobinary_char(output, data_bits_decomp, bytes_decomp); //.bc

  float* decompressed_data = readfrombinary_writetotxt_float(output, decomp, data_num);

  printf("done\n");
  free(bytes); 
  free(data);
  
  return 0;

}