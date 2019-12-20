/********************************************************************
 This program is for data compression used in MPI applications.
 ----------------------------------------------
 Email : huyao0107@gmail.com
 ---------------------------------------------------------------
********************************************************************/

#include <stdio.h>
#include "mpi.h"
#include "param.h"
#include "dataCompression.h"

float jacobi(int);
int initmax(int,int,int);
void initmt(int,int);
void initcomm(int,int,int);
void sendp(int,int,int);
void sendp1();
void sendp2();
void sendp3();

double fflop(int,int,int);
double mflops(int,double,double);

static float omega;
static int npe,id;

MPI_Datatype
myCompress_himeno(void* data, int count, int blklen, int stride, int starti, int startj, int startk)
{
  float real_value, predict_value1, predict_value2, predict_value3;
  int num = 0;

  for(int i=0; i<count; i++)
  {
    for(int j=0; j<blklen; j++)
    {
      real_value = data[starti][startj][startk+blklen]
    }
    if(stride == MKMAX)
    {
      startj++;
    }
    else if(stride == MJMAX*MKMAX)
    {
      starti++;
    }
  }

}
            
