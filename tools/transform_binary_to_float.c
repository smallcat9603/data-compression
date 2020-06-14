#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  char* float2 = argv[1]; //"00111100101011100010100100000000"; //"00111110011100000010110010000001";
  float float10 = strtofloat(float2);
  printf("%f \n", float10);

}