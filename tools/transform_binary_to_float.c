#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  char* float2 = "01000010110010000000000000000000";
  float float10 = strtofloat(float2);
  printf("%f \n", float10);

}