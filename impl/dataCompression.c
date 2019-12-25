/********************************************************************
 This program is for data compression used in MPI applications.
 ----------------------------------------------
 Email : huyao0107@gmail.com
 ---------------------------------------------------------------
********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include "param.h"
#include "dataCompression.h"

//myCompress
float calcCompressionRatio_himeno_ij_ik_jk(float data[MIMAX][MJMAX][MKMAX], int ijk, int v)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  float diff1, diff2, diff3, diff_min, selected_predict_value;
  int array_float_len = 0, array_char_len = 0;
  char compress_type;
  float compress_ratio;
  float* array_float = NULL; //(float*)malloc(sizeof(float));
  float* array_float_more = NULL;
  char* array_char = NULL; //(char*)malloc(sizeof(char));
  char* array_char_more = NULL;
  int* array_char_displacement = NULL;
  int* array_char_displacement_more = NULL;
  int A, B;

  if(ijk == 1) 
  {
    A = MJMAX;
    B = MKMAX;
  }
  else if(ijk == 2) 
  {
    A = MIMAX;
    B = MKMAX;
  }
  else if(ijk == 3)
  {
    A = MIMAX;
    B = MJMAX;
  } 
  for(int a=0; a<A; a++)
  {
    for(int b=0; b<B; b++)
    {
      if(ijk == 1) real_value = data[v][a][b]; // v is const
      else if(ijk == 2) real_value = data[a][v][b];
      else if(ijk == 3) real_value = data[a][b][v];

      if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
      {
        array_float_len++;
        array_float_more = (float*)realloc(array_float, sizeof(float)*array_float_len);
        if (array_float_more != NULL) 
        {
          array_float = array_float_more;
          array_float[array_float_len-1] = real_value;
        }
        else 
        {
          free(array_float);
          printf("Error (re)allocating memory");
          exit(1);
        }        
        
        if(before_value3 == -1) 
        {
          before_value3 = real_value; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = real_value;
        }
        else if(before_value1 == -1) 
        {
          before_value1 = real_value;
        }        
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
        compress_type = 'a';
        selected_predict_value = predict_value1;
        if(diff2<diff_min)
        {
          diff_min = diff2;
          compress_type = 'b';
          selected_predict_value = predict_value2;
        }
        if(diff3<diff_min)
        {
          diff_min = diff3;
          compress_type = 'c';
          selected_predict_value = predict_value3;
        }        

        before_value3 = before_value2;
        before_value2 = before_value1;
        if(diff_min<=absErrBound) 
        {
          array_char_len++;
          array_char_more = (char*)realloc(array_char, sizeof(char)*array_char_len);
          array_char_displacement_more = (int*)realloc(array_char_displacement, sizeof(int)*array_char_len);
          if (array_char_more != NULL && array_char_displacement_more != NULL) 
          {
            array_char = array_char_more;
            array_char[array_char_len-1] = compress_type;
            array_char_displacement = array_char_displacement_more;
            array_char_displacement[array_char_len-1] = array_float_len + array_char_len;
          }
          else 
          {
            free(array_char);
            free(array_char_displacement);
            printf("Error (re)allocating memory");
            exit(1);
          } 
          before_value1 = selected_predict_value;
        }
        else 
        {
          array_float_len++;
          array_float_more = (float*)realloc(array_float, sizeof(float)*array_float_len);
          if (array_float_more != NULL) 
          {
            array_float = array_float_more;
            array_float[array_float_len-1] = real_value;
          }
          else 
          {
            free(array_float);
            printf("Error (re)allocating memory");
            exit(1);
          }             
          before_value1 = real_value;
        }
      }
    }
  } 
  if(byte_or_bit == 1)
  {
    compress_ratio = (float)(array_char_len*sizeof(char)+array_float_len*sizeof(float))/((array_char_len+array_float_len)*sizeof(float));
  }
  else if(byte_or_bit == 2)
  {
    compress_ratio = (float)(array_char_len*2+array_float_len*sizeof(float)*8)/((array_char_len+array_float_len)*sizeof(float)*8);
  }
  return compress_ratio;
}

float calcCompressionRatio_himeno_sz(float data[MIMAX][MJMAX][MKMAX], int ijk, int v)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  float diff1, diff2, diff3, diff_min, predict_diff, selected_predict_value;
  char compress_type;
  float compress_ratio;
  int A, B;
  long origin_bits=0, compressed_bits=0;

  if(ijk == 1) 
  {
    A = MJMAX;
    B = MKMAX;
  }
  else if(ijk == 2) 
  {
    A = MIMAX;
    B = MKMAX;
  }
  else if(ijk == 3)
  {
    A = MIMAX;
    B = MJMAX;
  } 
  for(int a=0; a<A; a++)
  {
    for(int b=0; b<B; b++)
    {
      if(ijk == 1) real_value = data[v][a][b]; // v is const
      else if(ijk == 2) real_value = data[a][v][b];
      else if(ijk == 3) real_value = data[a][b][v];

      origin_bits += sizeof(float)*8;     

      if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
      { 
        compressed_bits += sizeof(float)*8; 
        if(before_value3 == -1) 
        {
          before_value3 = real_value; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = real_value;
        }
        else if(before_value1 == -1) 
        {
          before_value1 = real_value;
        }        
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
        compress_type = 'a';
        selected_predict_value = predict_value1;
        if(diff2<diff_min)
        {
          diff_min = diff2;
          compress_type = 'b';
          selected_predict_value = predict_value2;
        }
        if(diff3<diff_min)
        {
          diff_min = diff3;
          compress_type = 'c';
          selected_predict_value = predict_value3;
        }        

        before_value3 = before_value2;
        before_value2 = before_value1;
        if(diff_min<=absErrBound) 
        {
          if(byte_or_bit == 1)
          {
            compressed_bits += sizeof(char)*8; 
          }
          else if(byte_or_bit == 2)
          {
            compressed_bits += 2; 
          }
          before_value1 = selected_predict_value;
        }
        else 
        {
          float max, min;
          if(predict_value1 > predict_value2)
          {
            max = predict_value1;
            min = predict_value2;
          }
          else
          {
            max = predict_value2;
            min = predict_value1;
          }
          if(predict_value3 > max)
          {
            max = predict_value3;
          }
          else if(predict_value3 < min)
          {
            min = predict_value3;
          }
          
          predict_diff = max-min;

          char c[sizeof(float)*8];
          getFloatBin(predict_diff/2, c);
          int expo_value = 0;
          int mantissa_bits_within_error_bound;

          for(int i=1;i<9;i++) //1-9 exponential part of float (1-12 in the case of double)
          {
            if(c[i] != 0) 
            {
              expo_value += pow(2, 8-i);
            }  
          }
          expo_value -= 127;
          mantissa_bits_within_error_bound = absErrBound_binary + expo_value;
          if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
          {
            mantissa_bits_within_error_bound = 23;
          }
          else if(mantissa_bits_within_error_bound < 0)
          {
            mantissa_bits_within_error_bound = 0;
          }
          if(byte_or_bit == 1)
          {
            if(mantissa_bits_within_error_bound%8 != 0) compressed_bits += 1+8+(mantissa_bits_within_error_bound%8+1)*8;  
            else compressed_bits += 1+8+mantissa_bits_within_error_bound;
          }
          else if(byte_or_bit == 2)
          {
            compressed_bits += 1+8+mantissa_bits_within_error_bound;  
          }
          
          before_value1 = real_value;
        }
      }
    }
  } 
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_himeno_nolossy_performance(float data[MIMAX][MJMAX][MKMAX], int ijk, int v)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  float diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;
  int A, B;

  if(ijk == 1) 
  {
    A = MJMAX;
    B = MKMAX;
  }
  else if(ijk == 2) 
  {
    A = MIMAX;
    B = MKMAX;
  }
  else if(ijk == 3)
  {
    A = MIMAX;
    B = MJMAX;
  } 
  for(int a=0; a<A; a++)
  {
    for(int b=0; b<B; b++)
    {
      if(ijk == 1) real_value = data[v][a][b]; // v is const
      else if(ijk == 2) real_value = data[a][v][b];
      else if(ijk == 3) real_value = data[a][b][v];

      origin_bits += sizeof(float)*8;

      if(before_value4 == -1 || before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
      {
        compressed_bits += sizeof(float)*8;       
        
        if(before_value4 == -1) 
        {
          before_value4 = real_value; 
        }
        else if(before_value3 == -1) 
        {
          before_value3 = real_value; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = real_value;
        }
        else if(before_value1 == -1) 
        {
          before_value1 = real_value;
        }        
      }
      else
      {
        predict_value4 = 4*before_value1 - 6*before_value2 + 4*before_value3 - before_value4;

        diff4 = predict_value4-real_value;     

        before_value4 = before_value3;
        before_value3 = before_value2;
        before_value2 = before_value1;

        char c[sizeof(float)*8];
        getFloatBin(diff4, c);
        for(int i=1;i<sizeof(float)*8;i++)
        {
          if(c[i] != 0) 
          {
            if(byte_or_bit == 1)
            {
              if((sizeof(float)*8 - i + 3 + 1)%8 != 0) compressed_bits += ((sizeof(float)*8 - i + 3 + 1)%8+1)*8;  
              else compressed_bits += sizeof(float)*8 - i + 3 + 1;
            }
            else if(byte_or_bit == 2)
            {
              compressed_bits += sizeof(float)*8 - i + 3 + 1;
            }
            break;
          } 
        }
        before_value1 = real_value;
      }
    }
  }
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_himeno_nolossy_area(float data[MIMAX][MJMAX][MKMAX], int ijk, int v)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  float diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;
  int A, B;
  int cdb = 512, cdb_num = 1, occupied_bits = 0, indication = 5, le = 285, re = 222, re1 = 2, re2 = 4, re3 = 32, llrb = 2, ex = 1;

  if(ijk == 1) 
  {
    A = MJMAX;
    B = MKMAX;
  }
  else if(ijk == 2) 
  {
    A = MIMAX;
    B = MKMAX;
  }
  else if(ijk == 3)
  {
    A = MIMAX;
    B = MJMAX;
  } 
  for(int a=0; a<A; a++)
  {
    for(int b=0; b<B; b++)
    {
      if(ijk == 1) real_value = data[v][a][b]; // v is const
      else if(ijk == 2) real_value = data[a][v][b];
      else if(ijk == 3) real_value = data[a][b][v];

      origin_bits += sizeof(float)*8;

      if(before_value4 == -1 || before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
      {
        occupied_bits += re3+llrb+ex;       
        
        if(before_value4 == -1) 
        {
          before_value4 = real_value; 
        }
        else if(before_value3 == -1) 
        {
          before_value3 = real_value; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = real_value;
        }
        else if(before_value1 == -1) 
        {
          before_value1 = real_value;
        }        
      }
      else
      {
        predict_value4 = 4*before_value1 - 6*before_value2 + 4*before_value3 - before_value4;

        diff4 = predict_value4-real_value;     

        before_value4 = before_value3;
        before_value3 = before_value2;
        before_value2 = before_value1;

        char c[sizeof(float)*8];
        getFloatBin(diff4, c);
        for(int i=1;i<sizeof(float)*8;i++)
        {
          if(c[i] != 0) 
          {
            int nonzero = sizeof(float)*8 - i;
            int data_bits;
            if(nonzero <= re1)
            {
              data_bits = re1+llrb+ex;
            }
            else if(nonzero <= re2)
            {
              data_bits = re2+llrb+ex;
            }
            else if(nonzero <= re3)
            {
              data_bits = re3+llrb+ex;
            }
            
            if(occupied_bits + data_bits > cdb-indication)
            {
              cdb_num++;
              occupied_bits = data_bits;
            }
            else
            {
              occupied_bits += data_bits;
            }

            break;
          }
            
        }
        before_value1 = real_value;
      }
    }
  }
  compressed_bits = cdb_num*cdb;
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

void getFloatBin(float num,char bin[])
{
    int t = 1;//用来按位与操作
    int *f = (int*)(&num);//将float的解释成int，即float的地址转成int*
    for(int i=0;i<32;i++)
    {
    //从最高位开始按位与，如果为1，则bin[i]=1，如果为0，则bin[i]=0
    //这里没有将bin存成字符，而是数字1 0
        bin[i] = (*f)&(t<<31-i)?1:0;
    }
}

// MPI_Datatype
// myCompress_himeno(void* data, int count, int blklen, int stride, int starti, int startj, int startk)
// {
//   float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
//   float diff1, diff2, diff3, diff_min, selected_predict_value;
//   int array_float_len = 0, array_char_len = 0;
//   char compress_type;
//   float* array_float = NULL; //(float*)malloc(sizeof(float));
//   float* array_float_more = NULL;
//   char* array_char = NULL; //(char*)malloc(sizeof(char));
//   char* array_char_more = NULL;
//   int* array_char_displacement = NULL;
//   int* array_char_displacement_more = NULL;

//   for(int i=0; i<count; i++)
//   {
//     for(int j=0; j<blklen; j++)
//     {
//       real_value = data[starti][startj][startk+blklen];

//       if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
//       {
//         array_float_len++;
//         array_float_more = (float*)realloc(array_float, sizeof(float)*array_float_len);
//         if (array_float_more != NULL) 
//         {
//           array_float = array_float_more;
//           array_float[array_float_len-1] = real_value;
//         }
//         else 
//         {
//           free(array_float);
//           printf("Error (re)allocating memory");
//           exit(1);
//         }        
        
//         if(before_value3 == -1) 
//         {
//           before_value3 = real_value; 
//         }
//         else if(before_value2 == -1) 
//         {
//           before_value2 = real_value;
//         }
//         else if(before_value1 == -1) 
//         {
//           before_value1 = real_value;
//         }        
//       }
//       else
//       {
//         predict_value1 = before_value1;
//         predict_value2 = 2*before_value1 - before_value2;
//         predict_value3 = 3*before_value1 - 3*before_value2 + before_value3;

//         diff1 = fabs(predict_value1-real_value);
//         diff2 = fabs(predict_value2-real_value);
//         diff3 = fabs(predict_value3-real_value);

//         diff_min = diff1;
//         compress_type = 'a';
//         selected_predict_value = predict_value1;
//         if(diff2<diff_min)
//         {
//           diff_min = diff2;
//           compress_type = 'b';
//           selected_predict_value = predict_value2;
//         }
//         if(diff3<diff_min)
//         {
//           diff_min = diff3;
//           compress_type = 'c';
//           selected_predict_value = predict_value3;
//         }        

//         before_value3 = before_value2;
//         before_value2 = before_value1;
//         if(diff_min<=absErrBound) 
//         {
//           array_char_len++;
//           array_char_more = (char*)realloc(array_char, sizeof(char)*array_char_len);
//           array_char_displacement_more = (int*)realloc(array_char_displacement, sizeof(int)*array_char_len);
//           if (array_char_more != NULL && array_char_displacement_more != NULL) 
//           {
//             array_char = array_char_more;
//             array_char[array_char_len-1] = compress_type;
//             array_char_displacement = array_char_displacement_more;
//             array_char_displacement[array_char_len-1] = array_float_len + array_char_len;
//           }
//           else 
//           {
//             free(array_char);
//             free(array_char_displacement);
//             printf("Error (re)allocating memory");
//             exit(1);
//           } 
//           before_value1 = selected_predict_value;
//         }
//         else 
//         {
//           array_float_len++;
//           array_float_more = (float*)realloc(array_float, sizeof(float)*array_float_len);
//           if (array_float_more != NULL) 
//           {
//             array_float = array_float_more;
//             array_float[array_float_len-1] = real_value;
//           }
//           else 
//           {
//             free(array_float);
//             printf("Error (re)allocating memory");
//             exit(1);
//           }             
//           before_value1 = real_value;
//         }
//       }
//     }
//     if(stride == MKMAX)
//     {
//       startj++;
//     }
//     else if(stride == MJMAX*MKMAX)
//     {
//       starti++;
//     }
//   }
// }
            
