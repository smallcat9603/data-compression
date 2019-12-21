/********************************************************************
 This program is for data compression used in MPI applications.
 ----------------------------------------------
 Email : huyao0107@gmail.com
 ---------------------------------------------------------------
********************************************************************/

#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include "param.h"
#include "dataCompression.h"


MPI_Datatype myCompress_himeno(void*, int, int, int, int, int, int);

MPI_Datatype
myCompress_himeno(void* data, int count, int blklen, int stride, int starti, int startj, int startk)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  float diff1, diff2, diff3, diff_min, selected_predict_value;
  int array_float_len = 1, array_char_len = 1;
  float* array_float = (float*)malloc(sizeof(float)*array_float_len);
  char* array_char = (char*)malloc(sizeof(char)*array_char_len);

  for(int i=0; i<count; i++)
  {
    for(int j=0; j<blklen; j++)
    {
      real_value = data[starti][startj][startk+blklen];

      if(before_value3 == -1) 
      {
        before_value3 = real_value;
        array_float[array_float_len-1] = real_value;
        array_float_len++;
        float* array_float_old = array_float;
        array_float = (float*)realloc(sizeof(float)*array_float_len);
        if(array_float_old != array_float) free(array_float_old);
      }
      else if(before_value2 == -1) 
      {
        before_value2 = real_value;
      }
      else if(before_value1 == -1) 
      {
        before_value1 = real_value;
      }
      else
      {
        predict_value1 = before_value1;
        predict_value2 = 2*before_value1 - before_value2;
        predict_value3 = 3*before_value1 - 3*before_value2 + before_value3;

        diff1 = fabs(predict_value1-real_value);
        diff2 = fabs(predict_value2-real_value);
        diff3 = fabs(predict_value3-real_value);

        diff_min = diff1;
        selected_predict_value = predict_value1;
        if(diff2<diff_min)
        {
          diff_min = diff2;
          selected_predict_value = predict_value2;
        }
        if(diff3<diff_min)
        {
          diff_min = diff3;
          selected_predict_value = predict_value3;
        }        

        before_value3 = before_value2;
        before_value2 = before_value1;
        if(diff_min<=absErrBound) 
        {
          before_value1 = selected_predict_value;
        }
        else 
        {
          before_value1 = real_value;
        }
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
            
