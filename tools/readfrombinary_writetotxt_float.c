#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../impl/param.h"
#include "../impl/dataCompression.h"

int main(int argc, char** argv) {

  const char *binaryfile = "../spdp/num_plasma.sp";
  const char *txtfile = "../spdp/num_plasma.sp.txt";
  readfrombinary_writetotxt_float(binaryfile, txtfile, 4386200);

}