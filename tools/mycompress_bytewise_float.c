#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <time.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

struct vector
{
  float* p_data; //precise data
  char* c_data; //compressed data
};

int main(int argc, char** argv) {

  char input_txt[64], output_txt[128]; 

  if(argc < 2)
  {
    printf("Test case: ./mycompress_bytewise_float [input_txt]\n");
    printf("Example: ./mycompress_bytewise_float input.txt\n");
    exit(0);
  }
  
  sprintf(input_txt, "%s", argv[1]);
  printf("input_txt = %s\n", input_txt); 
  sprintf(output_txt, "%s.byte.txt", input_txt);

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

  float* array_float = NULL;
  char* array_char = NULL;
  int* array_char_displacement = NULL;

  double start_time_comp = clock();
  int array_float_len = myCompress(data, &array_float, &array_char, &array_char_displacement, num);
  double end_time_comp = clock();

  struct vector msg; 
  int num_p = array_float_len, num_c = num-array_float_len;
  msg.p_data = array_float;
  msg.c_data = array_char;

  double start_time_decomp = clock();
  float* decompressed_data = myDecompress(array_float, array_char, array_char_displacement, num);
  double end_time_decomp = clock();

  fp = fopen(output_txt, "w");
  for(int i=0; i<num; i++)
  {
    fprintf(fp, "%f\n", decompressed_data[i]);
  }  
  fclose(fp);
  printf("%sに保存しました。\n", output_txt);

  printf("absErrorBound: %f \n", absErrorBound); 

  float compress_ratio = (float)(num_c*sizeof(char)+num_p*sizeof(float))/((num_c+num_p)*sizeof(float));
  printf("Compression rate: %f \n", 1/compress_ratio); 
  printf("Compression time: %f sec \n", (end_time_comp-start_time_comp)/CLOCKS_PER_SEC); 
  printf("Decompression time: %f sec \n", (end_time_decomp-start_time_decomp)/CLOCKS_PER_SEC); 

  printf("done\n");

  free(data);
  
  return 0;

}