//
// Ping pong example with MPI_Send and MPI_Recv. Two processes ping pong a
// number back and forth, incrementing it until it reaches a given value.
//

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "param.h"
#include "dataCompression.h"

int main(int argc, char** argv) {
  // char* double2 = "0100000001011001000000000000000000000000000000000000000000000000";
  // char* float2 = "01000010110010000000000000000000";
  // //double2 = "0001010000000000100000000000000000111000011100100000001101111101"; //14 00 80 00 38 72 03 7D
  // double double10 = strtodbl(double2);
  // float float10 = strtofloat(float2);
  // printf("%f \n", double10);
  // printf("%f \n", float10);
  // double10 = 100;
  // float10 = 100;
  // char float_arr[32+1];
  // char double_arr[64+1];
  // floattostr(&float10, float_arr);
  // doubletostr(&double10, double_arr);
  // printf("%s \n", float_arr);
  // printf("%s \n", double_arr);  

  const char *file = "../fpc/test.trace";
  double *arr = readfrombinary_double(file, 100000);
  printf("%f \n", arr[0]);
  printf("%f \n", arr[1]);
  printf("%f \n", arr[2]);

}