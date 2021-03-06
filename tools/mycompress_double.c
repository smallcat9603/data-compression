#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  char input_txt[64], input_bi[64], output_bc[64], output_txt[64]; 

  if(argc < 2)
  {
    printf("Test case: mycompress_double [input_txt]\n");
    printf("Example: mycompress_double input.txt\n");
    exit(0);
  }
  
  sprintf(input_txt, "%s", argv[1]);
  printf("input_txt = %s\n", input_txt); 
  sprintf(input_bi, "%s.bi", input_txt);
  sprintf(output_bc, "%s.bc", input_txt);
  sprintf(output_txt, "%s.bc.txt", input_txt);

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

  unsigned char* data_bits_comp = NULL;
  int bytes_comp = 0; //total bytes of compressed data
  int pos = 8; //position of filled bit in last byte --> 87654321
  myCompress_bitwise_double(data_small, num, &data_bits_comp, &bytes_comp, &pos);
  writetobinary_char(output_bc, data_bits_comp, bytes_comp); //.bc

  int bytes_decomp = 0;
  unsigned char* data_bits_decomp = readfrombinary_char(output_bc, &bytes_decomp);

  double* decompressed_data = myDecompress_bitwise_double(data_bits_decomp, bytes_decomp, num);

  fp = fopen(output_txt, "w");
  for(int i=0; i<num; i++)
  {
    fprintf(fp, "%f\n", decompressed_data[i]+min);
  }  
  fclose(fp);
  printf("%sに保存しました。\n", output_txt);

  printf("absErrorBound: %f \n", absErrorBound); 

  float compress_ratio = (float)(bytes_comp*8)/(num*sizeof(double)*8);
  printf("Compression rate: %f \n", 1/compress_ratio); 

  printf("done\n");

  free(data);
  free(data_small);
  free(data_bits_comp);
  free(data_bits_decomp);
  free(decompressed_data);
  
  return 0;

}