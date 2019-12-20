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

MPI_Datatype myCompress_himeno(void*, int, int, int, int, int, int);

MPI_Datatype
myCompress_himeno(void* data, int count, int blklen, int stride, int starti, int startj, int startk)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  int num = 0;

  for(int i=0; i<count; i++)
  {
    for(int j=0; j<blklen; j++)
    {
      real_value = data[starti][startj][startk+blklen];

      if(before_value3 == -1) before_value3 = real_value;
      else if(before_value2 == -1) before_value2 = real_value;
      else if(before_value1 == -1) before_value1 = real_value;
      else
      {
        predict_value1 = before_value1;
        predict_value2 = 2*before_value1 - before_value2;
        predict_value3 = 3*before_value1 - 3*before_value2 + before_value3;

        

        before_value3 = before_value2;
        before_value2 = before_value1;
        before_value1 = 
      }
      
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
            
