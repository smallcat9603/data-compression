/********************************************************************
 This program is for data compression used in MPI applications.
 ----------------------------------------------
 Email : huyao0107@gmail.com
 ---------------------------------------------------------------
********************************************************************/

#include <stdio.h>
#include "mpi.h"
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

static float  p[MIMAX][MJMAX][MKMAX];
static float  a[4][MIMAX][MJMAX][MKMAX],
              b[3][MIMAX][MJMAX][MKMAX],
              c[3][MIMAX][MJMAX][MKMAX];
static float omega;
static int npe,id;

MPI_Datatype
myCompress(void* data, int count, int blklen, int stride, int start)
{
  float value;

  for(int i=0; i<count; i++)
  {
    for(int j=0; j<blklen; j++)
    {
      if(i+j>2)
      value = *(data+start+blklen+count*stride)
    }
  }

}
            
