#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  const char *binaryfile = "../fpc/num_plasma.trace";
  const char *txtfile = "../fpc/num_plasma.trace.txt";
  readfrombinary_writetotxt_double(binaryfile, txtfile, 4386200);

}