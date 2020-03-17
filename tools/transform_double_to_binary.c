#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  double double10 = atof(argv[1]);
  char double_arr[64+1];
  doubletostr(&double10, double_arr);
  printf("%s \n", double_arr);  

}