#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  char* double2 = argv[1]; //"0011111101001001101101000000000000000000000000000000000000000000";
  //double2 = "0001010000000000100000000000000000111000011100100000001101111101"; //14 00 80 00 38 72 03 7D
  double double10 = strtodbl(double2);
  printf("%.10lf \n", double10);

}