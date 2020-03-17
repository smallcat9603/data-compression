#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  float float10 = atof(argv[1]);
  char float_arr[32+1];
  floattostr(&float10, float_arr);
  printf("%s \n", float_arr);

}