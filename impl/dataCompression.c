/********************************************************************
 This program is for data compression used in MPI applications.
 ----------------------------------------------
 Email : huyao0107@gmail.com
 ---------------------------------------------------------------
********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <assert.h>
#include <zlib.h>
#include <time.h>
#include "mpi.h"
#include "param.h"
#include "dataCompression.h"

double absErrBound = absErrorBound;
int absErrorBound_binary = -100;

//double
void MPI_Bcast_bitwise_crc(double *buffer, int count, int root, int rank, int procs, float* compress_ratio, double* gosa, int* resend)
{
  uint32_t crc = 0;
  uint32_t crc_check = 0;
  unsigned char crc_ok = 'y';
  unsigned char* crc_ok_recv = NULL;  
  srand((unsigned)time(NULL));  

  int data_bytes = 0;
  double min = 0;
  unsigned char* data_bits = NULL;
  
  if(rank == root)
  {
      // sz_comp_ratio += calcCompressionRatio_sz_double(m_a, size_a);
      // nolossy_performance += calcCompressionRatio_nolossy_performance_double(m_a, size_a);
      // nolossy_area += calcCompressionRatio_nolossy_area_double(m_a, size_a);

      //mycommpress
      double* small = NULL;
      min = toSmallDataset_double(buffer, &small, count);

      int data_pos = 8; //position of filled bit in last byte --> 87654321

      myCompress_bitwise_double(small, count, &data_bits, &data_bytes, &data_pos);	
      crc = do_crc32(data_bits, data_bytes);		
  }

  MPI_Bcast(&data_bytes, 1, MPI_INT, root, MPI_COMM_WORLD);
  MPI_Bcast(&min, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
  *compress_ratio += data_bytes*8.0/(count*sizeof(double)*8);

  if(rank != root)
  {
      data_bits = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes);
  }
  MPI_Bcast(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);
  MPI_Bcast(&crc, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);

  if(rank != root)
  {
      crc_check = do_crc32(data_bits, data_bytes);

      if(BER > 0)
      {
          double ber = BER;
          uint64_t to = 1/ber;
          uint64_t r = get_random_int(0, to);
          if(r < data_bytes * 8)
          {
              crc_check = 0;
          }
      }
      
      if (crc == crc_check)
      {
          // printf("CRC passed\n");
          crc_ok = 'y';
      }  
      else
      {
          // printf("CRC NOT passed\n");
          crc_ok = 'n';
      }            
  }
  else
  {
      crc_ok_recv = (unsigned char *)malloc(procs*1*sizeof(unsigned char));
  }

  MPI_Gather(&crc_ok, 1, MPI_UNSIGNED_CHAR, crc_ok_recv, 1, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);
  
  if(rank == root)
  {
      for(int i = 0; i < procs; i++)
      {
          if(i != root && crc_ok_recv[i] == 'n')
          {
              MPI_Send(data_bits, data_bytes, MPI_UNSIGNED_CHAR, i, i, MPI_COMM_WORLD);
              (*resend)++;
          }
      }
  }
  else if(crc_ok == 'n')
  {
      MPI_Recv(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  double* decompressed_data = myDecompress_bitwise_double(data_bits, data_bytes, count);
  for(int i=0; i<count; i++)
  {
      if(rank == root)
      {
          *gosa += fabs(decompressed_data[i] + min - buffer[i]);
      }
      else
      {
          buffer[i] = decompressed_data[i] + min;
      }
  }

  // *gosa = *gosa/count;

  free(data_bits);
}

double* myDecompress_bitwise_double_mask(unsigned char* data_bits, int bytes, int num, int type, char mask[1+11+8])
{
  int offset_bits = 0;
  char* bits = NULL;
  char* bits_more = NULL;
  int bits_num = 0;
  int min_shift = 0;
  bool pending = false;

  double before_value1=-1, before_value2=-1, before_value3=-1;
  double* decompressed = (double*) malloc(sizeof(double)*num);
  int decompressed_num = 0;

  for(int i=0; i<bytes; i++)
  {
    for (int j=7; j>=min_shift; j--) //each bit of byte
    {
      int bit = (data_bits[i] >> j) & 1;

      if(offset_bits == 0) //start bit
      {
        if(bits_num == 0) //not start bit of mantissa
        {
          if(bit == 0)
          {
            pending = true;
            offset_bits = 1+type;
            //offset_bits = 1+8;          
          }
          else if(bit == 1)
          {
            offset_bits = 3; //100, 101, 110, 111            
          }          
        }
        else 
        {
          if(pending == true)
          {
            pending = false;
            bool masked = true;
            for (int n = 1; n < type+1; n++)
            {
              if(bits[n] != '1')
              {
                masked = false;
                break;
              }
            }
            if(masked == true)
            {
              offset_bits = 1;
            }
            else
            {
              offset_bits = 11 - type;
            }          
          }
          else //start bit of mantissa
          {
            int expo_value = 0;
            for(int n=1; n<12; n++)
            {
              if(bits_num == 1+11)
              {
                expo_value += (bits[n]-'0')*pow(2,11-n);
              }
              else if(bits_num == 1+type+1)
              {
                expo_value += (mask[n]-'0')*pow(2,11-n);
              }
              else
              {
                printf("bits_num error");
                exit(1);
              }
            }
            expo_value -= 1023;           

            if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

            int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
            if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
            {
              mantissa_bits_within_error_bound = 52;
            }
            else if(mantissa_bits_within_error_bound < 0)
            {
              mantissa_bits_within_error_bound = 0;
            }

            offset_bits = mantissa_bits_within_error_bound;

            if(offset_bits > 0) //has mantissa bits
            {
              if(bits_num == 1 + type + 1) //0 11 0/1
              {
                if(bits[bits_num-1] == '0')
                {
                  offset_bits -= 8; //8 --> 0
                }
                // else if(bits[bits_num-1] == '1')
                // {
                //   offset_bits -= 4; //8 --> 4
                // }
              }
            }
            else //no mantissa bit
            {
              decompressed_num++;
              decompressed[decompressed_num-1] = decompress_bitwise_double_mask(bits, bits_num, before_value1, before_value2, before_value3, type, mask);

              if(before_value3 == -1) 
              {
                before_value3 = decompressed[decompressed_num-1]; 
              }
              else if(before_value2 == -1) 
              {
                before_value2 = decompressed[decompressed_num-1];
              }
              else if(before_value1 == -1) 
              {
                before_value1 = decompressed[decompressed_num-1];
              }
              else
              {
                before_value3 = before_value2;
                before_value2 = before_value1;
                before_value1 = decompressed[decompressed_num-1];
              }

              bits = NULL;
              bits_num = 0;
              pending = false;

              if(bit == 0)
              {
                pending = true;
                offset_bits = 1+type;              
                //offset_bits = 1+8;            
              }
              else if(bit == 1)
              {
                offset_bits = 3;             
              }                           
            }                        
          }          
        }      
      }
      bits_num++;
      bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
      if (bits_more != NULL) 
      {
        bits = bits_more;
        bits[bits_num-1] = bit + '0';
      }
      else 
      {
        free(bits);
        printf("Error (re)allocating memory\n");
        exit(1);
      }   

      offset_bits--;
      if(offset_bits == 0 && bits_num != 1+11 && bits_num != 1+type+1 && pending == false)
      {
        decompressed_num++;
        decompressed[decompressed_num-1] = decompress_bitwise_double_mask(bits, bits_num, before_value1, before_value2, before_value3, type, mask);
        // printf("%f ", decompressed[decompressed_num-1]);
        
        if(before_value3 == -1) 
        {
          before_value3 = decompressed[decompressed_num-1]; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = decompressed[decompressed_num-1];
        }
        else if(before_value1 == -1) 
        {
          before_value1 = decompressed[decompressed_num-1];
        }
        else
        {
          before_value3 = before_value2;
          before_value2 = before_value1;
          before_value1 = decompressed[decompressed_num-1];
        }

        bits = NULL;
        bits_num = 0;
        pending = false;
      }
    }       
  }
  return decompressed;
}

double decompress_bitwise_double_mask(char* bits, int bits_num, double before_value1, double before_value2, double before_value3, int type, char mask[1+11+8])
{
  if(bits_num == 3 && bits[0] == '1')
  {
    if(bits[1] == '0' && bits[2] == '0')
    {
      return 0.0;
    }
    else if(bits[1] == '0' && bits[2] == '1')
    {
      return before_value1;
    }
    else if(bits[1] == '1' && bits[2] == '0')
    {
      return 2*before_value1 - before_value2;
    }
    else if(bits[1] == '1' && bits[2] == '1')
    {
      return 3*before_value1 - 3*before_value2 + before_value3;
    }
  }
  else
  {
    if(bits_num == sizeof(double)*8)
    {
      return strtodbl(bits);
    }
    else
    {
      bool masked = true;
      char* bits64 = NULL;
      for(int n=1; n<type+1; n++)
      {
        if(bits[n] != '1')
        {
          masked = false;
          break;
        }
      }
      if(masked == true)
      {
        // bits32 = (char*)realloc(mask, sizeof(float)*8);
        bits64 = malloc(sizeof(double)*8);
        for(int i=0; i<1+11+8; i++)
        {
          bits64[i] = mask[i];
        }
        if(bits[type+1] == '0')
        {
          for(int i=20, j=1+type+1; j<bits_num; i++, j++)
          {
            bits64[i] = bits[j];
          }

          bits64[20+bits_num-(1+type+1)] = '1';
          if(20+bits_num-(1+type+1)+1 < sizeof(double)*8)     
          {
            for(int i=20+bits_num-(1+type+1)+1; i< sizeof(double)*8; i++)
            {
              bits64[i] = '0';
            }
          }          
        }
        else if(bits[type+1] == '1')
        {
          for(int i=12, j=1+type+1; j<bits_num; i++, j++)
          {
            bits64[i] = bits[j];
          }          
   
          bits64[12+bits_num-(1+type+1)] = '1';
          if(12+bits_num-(1+type+1)+1 < sizeof(double)*8)     
          {
            for(int i=12+bits_num-(1+type+1)+1; i< sizeof(double)*8; i++)
            {
              bits64[i] = '0';
            }
          }                                          
        }
      }
      else
      {
        bits64 = (char*)realloc(bits, sizeof(double)*8);
        bits64[bits_num] = '1';
        if(bits_num+1 < sizeof(double)*8)     
        {
          for(int i=bits_num+1; i< sizeof(double)*8; i++)
          {
            bits64[i] = '0';
          }
        }
      }
      return strtodbl(bits64);   
    }
  }
}

void compress_bitwise_double_mask(double real_value, unsigned char** data_bits, int* bytes, int* pos, int type, char mask[1+11+8])
{
  double double10 = real_value;
  char double_arr[64+1];
  doubletostr(&double10, double_arr);

  int expo_value = 0;
  for(int i=1; i<12; i++)
  {
    expo_value += (double_arr[i]-'0')*pow(2,11-i);
  }
  expo_value -= 1023;  

  if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

  int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;

  if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
  {
    mantissa_bits_within_error_bound = 52;
  }
  else if(mantissa_bits_within_error_bound < 0)
  {
    mantissa_bits_within_error_bound = 0;
  }

  int bits_after_compress = 1+11+mantissa_bits_within_error_bound;  

  //bitmask
  bool masked = true;
  for(int i=0; i<12; i++)
  {
    if(double_arr[i] != mask[i])
    {
      masked = false;
      break;
    }
  }

  if(masked == true)
  {
    int error = 0;
    int index = -1;

    for(int i=12; i<20; i++)
    {
      if(double_arr[i] != mask[i])
      {
        error++;
        break;
      }  
    }  

    if(error == 0)
    {
      //type=1 --> 01, type=2 --> 011, type=3 --> 0111, ...
      add_bit_to_bytes(data_bits, bytes, pos, 0);
      for(int i=0; i<type; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
      }
      add_bit_to_bytes(data_bits, bytes, pos, 0);   
      for(int i=20; i<bits_after_compress; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, double_arr[i]-'0');
      }       
    }
    else if(error == 1)
    {
      //type=1 --> 01, type=2 --> 011, type=3 --> 0111, ...
      add_bit_to_bytes(data_bits, bytes, pos, 0);
      for(int i=0; i<type; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
      }
      add_bit_to_bytes(data_bits, bytes, pos, 1);
       
      for(int i=12; i<bits_after_compress; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, double_arr[i]-'0');
      }    
    }
    else
    {
      printf("error error");
      exit(0);      
    }  
  }
  else
  {
    for(int i=0; i<bits_after_compress; i++)
    {
      add_bit_to_bytes(data_bits, bytes, pos, double_arr[i]-'0');
    }
  }
}

void myCompress_bitwise_double_mask(double data[], int num, unsigned char** data_bits, int* bytes, int* pos, int type, char mask[1+11+8])
{
  double real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  double diff1, diff2, diff3, diff_min, selected_predict_value;
  char compress_type;
  int a=0, b=0, c=0, d=0;

  for(int n=0; n<num; n++)
  {
    real_value = data[n];

    if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    {
      //if(real_value == 0)
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
      }
      else
      {
        compress_bitwise_double_mask(real_value, data_bits, bytes, pos, type, mask);
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
      compress_type = 'a'; //101
      selected_predict_value = predict_value1;
      if(diff2<diff_min)
      {
        diff_min = diff2;
        compress_type = 'b'; //110
        selected_predict_value = predict_value2;
      }
      if(diff3<diff_min)
      {
        diff_min = diff3;
        compress_type = 'c'; //111
        selected_predict_value = predict_value3;
      }        

      before_value3 = before_value2;
      before_value2 = before_value1;
      before_value1 = real_value;
      
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
      }
      else if(diff_min<=absErrorBound) 
      {
        if(compress_type == 'a')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);
          add_bit_to_bytes(data_bits, bytes, pos, 1);   
          a++;     
        }
        else if(compress_type == 'b')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);  
          b++;
        }
        else if(compress_type == 'c')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);  
          c++;
        }
        else
        {
          printf("Error compress_type\n");
          exit(1);
        }
      }
      else 
      {
        compress_bitwise_double_mask(real_value, data_bits, bytes, pos, type, mask);            
      }
    }
  }   
}

float* myDecompress_bitwise_mask(unsigned char* data_bits, int bytes, int num, int type, char mask[1+8+8])
{
  int offset_bits = 0;
  char* bits = NULL;
  char* bits_more = NULL;
  int bits_num = 0;
  int min_shift = 0;
  bool pending = false;

  float before_value1=-1, before_value2=-1, before_value3=-1;
  float* decompressed = (float*) malloc(sizeof(float)*num);
  int decompressed_num = 0;

  for(int i=0; i<bytes; i++)
  {
    for (int j=7; j>=min_shift; j--) //each bit of byte
    {
      int bit = (data_bits[i] >> j) & 1;

      if(offset_bits == 0) //start bit
      {
        if(bits_num == 0) //not start bit of mantissa
        {
          if(bit == 0)
          {
            pending = true;
            offset_bits = 1+type;
            //offset_bits = 1+8;          
          }
          else if(bit == 1)
          {
            offset_bits = 3; //100, 101, 110, 111            
          }          
        }
        else 
        {
          if(pending == true)
          {
            pending = false;
            bool masked = true;
            for (int n = 1; n < type+1; n++)
            {
              if(bits[n] != '1')
              {
                masked = false;
                break;
              }
            }
            if(masked == true)
            {
              offset_bits = 1;
            }
            else
            {
              offset_bits = 8 - type;
            }          
          }
          else //start bit of mantissa
          {
            int expo_value = 0;
            for(int n=1; n<9; n++)
            {
              if(bits_num == 1+8)
              {
                expo_value += (bits[n]-'0')*pow(2,8-n);
              }
              else if(bits_num == 1+type+1)
              {
                expo_value += (mask[n]-'0')*pow(2,8-n);
              }
              else
              {
                printf("bits_num error");
                exit(1);
              }
            }
            expo_value -= 127;           

            if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

            int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
            if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
            {
              mantissa_bits_within_error_bound = 23;
            }
            else if(mantissa_bits_within_error_bound < 0)
            {
              mantissa_bits_within_error_bound = 0;
            }

            offset_bits = mantissa_bits_within_error_bound;

            if(offset_bits > 0) //has mantissa bits
            {
              if(bits_num == 1 + type + 1) //0 11 0/1
              {
                if(bits[bits_num-1] == '0')
                {
                  offset_bits -= 8; //8 --> 0
                }
                // else if(bits[bits_num-1] == '1')
                // {
                //   offset_bits -= 4; //8 --> 4
                // }
              }
            }
            else //no mantissa bit
            {
              decompressed_num++;
              decompressed[decompressed_num-1] = decompress_bitwise_float_mask(bits, bits_num, before_value1, before_value2, before_value3, type, mask);

              if(before_value3 == -1) 
              {
                before_value3 = decompressed[decompressed_num-1]; 
              }
              else if(before_value2 == -1) 
              {
                before_value2 = decompressed[decompressed_num-1];
              }
              else if(before_value1 == -1) 
              {
                before_value1 = decompressed[decompressed_num-1];
              }
              else
              {
                before_value3 = before_value2;
                before_value2 = before_value1;
                before_value1 = decompressed[decompressed_num-1];
              }

              bits = NULL;
              bits_num = 0;
              pending = false;

              if(bit == 0)
              {
                pending = true;
                offset_bits = 1+type;              
                //offset_bits = 1+8;            
              }
              else if(bit == 1)
              {
                offset_bits = 3;             
              }                           
            }                        
          }          
        }      
      }
      bits_num++;
      bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
      if (bits_more != NULL) 
      {
        bits = bits_more;
        bits[bits_num-1] = bit + '0';
      }
      else 
      {
        free(bits);
        printf("Error (re)allocating memory\n");
        exit(1);
      }   

      offset_bits--;
      if(offset_bits == 0 && bits_num != 1+8 && bits_num != 1+type+1 && pending == false)
      {
        decompressed_num++;
        decompressed[decompressed_num-1] = decompress_bitwise_float_mask(bits, bits_num, before_value1, before_value2, before_value3, type, mask);
        // printf("%f ", decompressed[decompressed_num-1]);
        
        if(before_value3 == -1) 
        {
          before_value3 = decompressed[decompressed_num-1]; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = decompressed[decompressed_num-1];
        }
        else if(before_value1 == -1) 
        {
          before_value1 = decompressed[decompressed_num-1];
        }
        else
        {
          before_value3 = before_value2;
          before_value2 = before_value1;
          before_value1 = decompressed[decompressed_num-1];
        }

        bits = NULL;
        bits_num = 0;
        pending = false;
      }
    }       
  }
  return decompressed;
}

float decompress_bitwise_float_mask(char* bits, int bits_num, float before_value1, float before_value2, float before_value3, int type, char mask[1+8+8])
{
  if(bits_num == 3 && bits[0] == '1')
  {
    if(bits[1] == '0' && bits[2] == '0')
    {
      return 0.0;
    }
    else if(bits[1] == '0' && bits[2] == '1')
    {
      return before_value1;
    }
    else if(bits[1] == '1' && bits[2] == '0')
    {
      return 2*before_value1 - before_value2;
    }
    else if(bits[1] == '1' && bits[2] == '1')
    {
      return 3*before_value1 - 3*before_value2 + before_value3;
    }
  }
  else
  {
    if(bits_num == sizeof(float)*8)
    {
      return strtofloat(bits);
    }
    else
    {
      bool masked = true;
      char* bits32 = NULL;
      for(int n=1; n<type+1; n++)
      {
        if(bits[n] != '1')
        {
          masked = false;
          break;
        }
      }
      if(masked == true)
      {
        // bits32 = (char*)realloc(mask, sizeof(float)*8);
        bits32 = malloc(sizeof(float)*8);
        for(int i=0; i<1+8+8; i++)
        {
          bits32[i] = mask[i];
        }
        if(bits[type+1] == '0')
        {
          for(int i=17, j=1+type+1; j<bits_num; i++, j++)
          {
            bits32[i] = bits[j];
          }

          bits32[17+bits_num-(1+type+1)] = '1';
          if(17+bits_num-(1+type+1)+1 < sizeof(float)*8)     
          {
            for(int i=17+bits_num-(1+type+1)+1; i< sizeof(float)*8; i++)
            {
              bits32[i] = '0';
            }
          }          
        }
        else if(bits[type+1] == '1')
        {
          // if(bits[type+2] == '0' && bits[type+3] == '0')
          // {
          //   bits32[9] = bits[type+4];
          //   bits32[10] = bits[type+5];
          // }
          // else if(bits[type+2] == '0' && bits[type+3] == '1')
          // {
          //   bits32[11] = bits[type+4];
          //   bits32[12] = bits[type+5];
          // }   
          // else if(bits[type+2] == '1' && bits[type+3] == '0')
          // {
          //   bits32[13] = bits[type+4];
          //   bits32[14] = bits[type+5];
          // }     
          // else if(bits[type+2] == '1' && bits[type+3] == '1')
          // {
          //   bits32[15] = bits[type+4];
          //   bits32[16] = bits[type+5];
          // }
          // for(int i=17, j=1+type+1+4; j<bits_num; i++, j++)
          // {
          //   bits32[i] = bits[j];
          // }
          for(int i=9, j=1+type+1; j<bits_num; i++, j++)
          {
            bits32[i] = bits[j];
          }          

          // bits32[17+bits_num-(1+type+1+4)] = '1';
          // if(17+bits_num-(1+type+1+4)+1 < sizeof(float)*8)     
          // {
          //   for(int i=17+bits_num-(1+type+1+4)+1; i< sizeof(float)*8; i++)
          //   {
          //     bits32[i] = '0';
          //   }
          // }    
          bits32[9+bits_num-(1+type+1)] = '1';
          if(9+bits_num-(1+type+1)+1 < sizeof(float)*8)     
          {
            for(int i=9+bits_num-(1+type+1)+1; i< sizeof(float)*8; i++)
            {
              bits32[i] = '0';
            }
          }                                          
        }
      }
      else
      {
        bits32 = (char*)realloc(bits, sizeof(float)*8);
        bits32[bits_num] = '1';
        if(bits_num+1 < sizeof(float)*8)     
        {
          for(int i=bits_num+1; i< sizeof(float)*8; i++)
          {
            bits32[i] = '0';
          }
        }
      }
      return strtofloat(bits32);   
    }
  }
}

//bitmask-based bitwise myCompress for ping-pong & himeno (float)
void myCompress_bitwise_mask(float data[], int num, unsigned char** data_bits, int* bytes, int* pos, int type, char mask[1+8+8])
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  float diff1, diff2, diff3, diff_min, selected_predict_value;
  char compress_type;
  int a=0, b=0, c=0, d=0;

  for(int n=0; n<num; n++)
  {
    real_value = data[n];

    if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    {
      //if(real_value == 0)
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
      }
      else
      {
        compress_bitwise_float_mask(real_value, data_bits, bytes, pos, type, mask);
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
      compress_type = 'a'; //101
      selected_predict_value = predict_value1;
      if(diff2<diff_min)
      {
        diff_min = diff2;
        compress_type = 'b'; //110
        selected_predict_value = predict_value2;
      }
      if(diff3<diff_min)
      {
        diff_min = diff3;
        compress_type = 'c'; //111
        selected_predict_value = predict_value3;
      }        

      before_value3 = before_value2;
      before_value2 = before_value1;
      before_value1 = real_value;
      
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
      }
      else if(diff_min<=absErrorBound) 
      {
        if(compress_type == 'a')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);
          add_bit_to_bytes(data_bits, bytes, pos, 1);   
          a++;     
        }
        else if(compress_type == 'b')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);  
          b++;
        }
        else if(compress_type == 'c')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);  
          c++;
        }
        else
        {
          printf("Error compress_type\n");
          exit(1);
        }
      }
      else 
      {
        compress_bitwise_float_mask(real_value, data_bits, bytes, pos, type, mask);            
      }
    }
  }  
}

void compress_bitwise_float_mask(float real_value, unsigned char** data_bits, int* bytes, int* pos, int type, char mask[1+8+8])
{
  float float10 = real_value;
  char float_arr[32+1];
  floattostr(&float10, float_arr);

  int expo_value = 0;
  for(int i=1; i<9; i++)
  {
    expo_value += (float_arr[i]-'0')*pow(2,8-i);
  }
  expo_value -= 127;  

  if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

  int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;

  if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
  {
    mantissa_bits_within_error_bound = 23;
  }
  else if(mantissa_bits_within_error_bound < 0)
  {
    mantissa_bits_within_error_bound = 0;
  }

  int bits_after_compress = 1+8+mantissa_bits_within_error_bound;  

  //bitmask
  bool masked = true;
  for(int i=0; i<9; i++)
  {
    if(float_arr[i] != mask[i])
    {
      masked = false;
      break;
    }
  }

  if(masked == true)
  {
    int error = 0;
    int index = -1;
    // for(int i=9; i<17; i=i+2)
    // {
    //   if(float_arr[i] != mask[i] || float_arr[i+1] != mask[i+1])
    //   {
    //     error++;
    //     index = i; //9, 11, 13, 15
    //     if(error > 1)
    //     {
    //       index = -1;
    //       break;
    //     }
    //   }
    // }
    for(int i=9; i<17; i++)
    {
      if(float_arr[i] != mask[i])
      {
        error++;
        break;
      }  
    }  

    if(error == 0)
    {
      //type=1 --> 01, type=2 --> 011, type=3 --> 0111, ...
      add_bit_to_bytes(data_bits, bytes, pos, 0);
      for(int i=0; i<type; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
      }
      add_bit_to_bytes(data_bits, bytes, pos, 0);   
      for(int i=17; i<bits_after_compress; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, float_arr[i]-'0');
      }       
    }
    else if(error == 1)
    {
      //type=1 --> 01, type=2 --> 011, type=3 --> 0111, ...
      add_bit_to_bytes(data_bits, bytes, pos, 0);
      for(int i=0; i<type; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
      }
      add_bit_to_bytes(data_bits, bytes, pos, 1);
      // if(index == 9)
      // {
      //   add_bit_to_bytes(data_bits, bytes, pos, 0);
      //   add_bit_to_bytes(data_bits, bytes, pos, 0);
      // }
      // else if(index == 11)
      // {
      //   add_bit_to_bytes(data_bits, bytes, pos, 0);
      //   add_bit_to_bytes(data_bits, bytes, pos, 1);
      // }
      // else if(index == 13)
      // {
      //   add_bit_to_bytes(data_bits, bytes, pos, 1);
      //   add_bit_to_bytes(data_bits, bytes, pos, 0);          
      // }
      // else if(index == 15)
      // {
      //   add_bit_to_bytes(data_bits, bytes, pos, 1);
      //   add_bit_to_bytes(data_bits, bytes, pos, 1);             
      // }
      // else
      // {
      //   printf("index error");
      //   exit(0);
      // }
      // add_bit_to_bytes(data_bits, bytes, pos, float_arr[index]-'0');
      // add_bit_to_bytes(data_bits, bytes, pos, float_arr[index+1]-'0');          
      // for(int i=17; i<bits_after_compress; i++)
      // {
      //   add_bit_to_bytes(data_bits, bytes, pos, float_arr[i]-'0');
      // }        
      for(int i=9; i<bits_after_compress; i++)
      {
        add_bit_to_bytes(data_bits, bytes, pos, float_arr[i]-'0');
      }    
    }
    else
    {
      // for(int i=0; i<bits_after_compress; i++)
      // {
      //   add_bit_to_bytes(data_bits, bytes, pos, float_arr[i]-'0');
      // }
      printf("error error");
      exit(0);      
    }  
  }
  else
  {
    for(int i=0; i<bits_after_compress; i++)
    {
      add_bit_to_bytes(data_bits, bytes, pos, float_arr[i]-'0');
    }
  }
}

double* myDecompress_bitwise_double_np(unsigned char* data_bits, int bytes, int num)
{
  int offset_bits = 0;
  char* bits = NULL;
  char* bits_more = NULL;
  int bits_num = 0;
  int min_shift = 0;

  double* decompressed = (double*) malloc(sizeof(double)*num);
  int decompressed_num = 0;

  for(int i=0; i<bytes; i++)
  {
    //if(i == bytes - 1 && pos != 8) min_shift = pos;

    for (int j=7; j>=min_shift; j--) //each bit of byte
    {
      int bit = (data_bits[i] >> j) & 1;

      //printf("%d(%d)", bit, bits_num);

      if(offset_bits == 0) //start bit
      {
        if(bits_num == 0) //not start bit of mantissa
        {
          if(bit == 0)
          {
            offset_bits = 1+11;
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            }             
          }
          else if(bit == 1)
          {
            printf("Error leading bit 1\n");
            exit(1);             
          }
        }
        else //start bit of mantissa
        {
          int expo_value = 0;
          for(int i=1; i<12; i++)
          {
            expo_value += (bits[i]-'0')*pow(2,11-i);
          }
          expo_value -= 1023;

          if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

          int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
          if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
          {
            mantissa_bits_within_error_bound = 52;
          }
          else if(mantissa_bits_within_error_bound < 0)
          {
            mantissa_bits_within_error_bound = 0;
          }

          offset_bits = mantissa_bits_within_error_bound;

          if(offset_bits > 0) //has mantissa bits
          {
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            } 
          }
          else //no mantissa bit
          {
            decompressed_num++;
            decompressed[decompressed_num-1] = decompress_bitwise_double_np(bits, bits_num);
            // printf("%f ", decompressed[decompressed_num-1]);

            bits = NULL;
            bits_num = 0;

            if(bit == 0)
            {
              offset_bits = 1+11;
              bits_num++;
              bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
              if (bits_more != NULL) 
              {
                bits = bits_more;
                bits[bits_num-1] = bit + '0';
              }
              else 
              {
                free(bits);
                printf("Error (re)allocating memory\n");
                exit(1);
              }             
            }
            else if(bit == 1)
            {
              printf("Error leading bit 1\n");
              exit(1);            
            }              
          }
        }
      }
      else
      {
        bits_num++;
        bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
        if (bits_more != NULL) 
        {
          bits = bits_more;
          bits[bits_num-1] = bit + '0';
        }
        else 
        {
          free(bits);
          printf("Error (re)allocating memory\n");
          exit(1);
        }        
      }
      offset_bits--;
      if(offset_bits == 0 && bits_num != 1+11)
      {
        decompressed_num++;
        decompressed[decompressed_num-1] = decompress_bitwise_double_np(bits, bits_num);
        // printf("%f ", decompressed[decompressed_num-1]);

        bits = NULL;
        bits_num = 0;
      }
    }       
  }
  return decompressed;
}

double decompress_bitwise_double_np(char* bits, int bits_num)
{
  if(bits_num == sizeof(double)*8)
  {
    return strtodbl(bits);
  }
  else
  {
    char* bits64 = (char*)realloc(bits, sizeof(double)*8);
    bits64[bits_num] = '1';
    if(bits_num+1 < sizeof(double)*8)     
    {
      for(int i=bits_num+1; i< sizeof(double)*8; i++)
      {
        bits64[i] = '0';
      }
    }
    return strtodbl(bits64); 
  }
}

float* myDecompress_bitwise_np(unsigned char* data_bits, int bytes, int num)
{
  int offset_bits = 0;
  char* bits = NULL;
  char* bits_more = NULL;
  int bits_num = 0;
  int min_shift = 0;

  float* decompressed = (float*) malloc(sizeof(float)*num);
  int decompressed_num = 0;

  for(int i=0; i<bytes; i++)
  {
    //if(i == bytes - 1 && pos != 8) min_shift = pos;

    for (int j=7; j>=min_shift; j--) //each bit of byte
    {
      int bit = (data_bits[i] >> j) & 1;

      //printf("%d(%d)", bit, bits_num);

      if(offset_bits == 0) //start bit
      {
        if(bits_num == 0) //not start bit of mantissa
        {
          if(bit == 0)
          {
            offset_bits = 1+8;
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            }             
          }
          else if(bit == 1)
          {
            printf("Error leading bit 1\n");
            exit(1);            
          }
        }
        else //start bit of mantissa
        {
          int expo_value = 0;
          for(int i=1; i<9; i++)
          {
            expo_value += (bits[i]-'0')*pow(2,8-i);
          }
          expo_value -= 127;      

          if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

          int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
          if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
          {
            mantissa_bits_within_error_bound = 23;
          }
          else if(mantissa_bits_within_error_bound < 0)
          {
            mantissa_bits_within_error_bound = 0;
          }

          offset_bits = mantissa_bits_within_error_bound;

          if(offset_bits > 0) //has mantissa bits
          {
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            } 
          }
          else //no mantissa bit
          {
            decompressed_num++;
            decompressed[decompressed_num-1] = decompress_bitwise_float_np(bits, bits_num);
            // printf("%f ", decompressed[decompressed_num-1]);

            bits = NULL;
            bits_num = 0;

            if(bit == 0)
            {
              offset_bits = 1+8;
              bits_num++;
              bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
              if (bits_more != NULL) 
              {
                bits = bits_more;
                bits[bits_num-1] = bit + '0';
              }
              else 
              {
                free(bits);
                printf("Error (re)allocating memory\n");
                exit(1);
              }             
            }
            else if(bit == 1)
            {
              printf("Error leading bit 1\n");
              exit(1);            
            }              
          }
        }
      }
      else
      {
        bits_num++;
        bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
        if (bits_more != NULL) 
        {
          bits = bits_more;
          bits[bits_num-1] = bit + '0';
        }
        else 
        {
          free(bits);
          printf("Error (re)allocating memory\n");
          exit(1);
        }        
      }
      offset_bits--;
      if(offset_bits == 0 && bits_num != 1+8)
      {
        decompressed_num++;
        decompressed[decompressed_num-1] = decompress_bitwise_float_np(bits, bits_num);
        // printf("%f ", decompressed[decompressed_num-1]);

        bits = NULL;
        bits_num = 0;
      }
    }       
  }
  return decompressed;
}

float decompress_bitwise_float_np(char* bits, int bits_num)
{
  if(bits_num == sizeof(float)*8)
  {
    return strtofloat(bits);
  }
  else
  {
    char* bits32 = (char*)realloc(bits, sizeof(float)*8);
    bits32[bits_num] = '1';
    if(bits_num+1 < sizeof(float)*8)     
    {
      for(int i=bits_num+1; i< sizeof(float)*8; i++)
      {
        bits32[i] = '0';
      }
    }
    return strtofloat(bits32); 
  }
}

//bitwise myCompress no prediction for k-means (double)
void myCompress_bitwise_double_np(double data[], int num, unsigned char** data_bits, int* bytes, int* pos)
{
  double real_value;

  for(int n=0; n<num; n++)
  {
    real_value = data[n];
    compress_bitwise_double(real_value, data_bits, bytes, pos);
  }
}

//bitwise myCompress no prediction for ping-pong & himeno (float)
void myCompress_bitwise_np(float data[], int num, unsigned char** data_bits, int* bytes, int* pos)
{
  float real_value;

  for(int n=0; n<num; n++)
  {
    real_value = data[n];
    compress_bitwise_float(real_value, data_bits, bytes, pos);
  }
}

double* myDecompress_bitwise_double(unsigned char* data_bits, int bytes, int num)
{
  int offset_bits = 0;
  char* bits = NULL;
  char* bits_more = NULL;
  int bits_num = 0;
  int min_shift = 0;

  double before_value1=-1, before_value2=-1, before_value3=-1;
  double* decompressed = (double*) malloc(sizeof(double)*num);
  int decompressed_num = 0;

  for(int i=0; i<bytes; i++)
  {
    //if(i == bytes - 1 && pos != 8) min_shift = pos;

    for (int j=7; j>=min_shift; j--) //each bit of byte
    {
      int bit = (data_bits[i] >> j) & 1;

      //printf("%d(%d)", bit, bits_num);

      if(offset_bits == 0) //start bit
      {
        if(bits_num == 0) //not start bit of mantissa
        {
          if(bit == 0)
          {
            offset_bits = 1+11;
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            }             
          }
          else if(bit == 1)
          {
            offset_bits = 3; //100, 101, 110, 111
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            }             
          }
        }
        else //start bit of mantissa
        {
          int expo_value = 0;
          for(int i=1; i<12; i++)
          {
            expo_value += (bits[i]-'0')*pow(2,11-i);
          }
          expo_value -= 1023;            

          if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

          int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
          if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
          {
            mantissa_bits_within_error_bound = 52;
          }
          else if(mantissa_bits_within_error_bound < 0)
          {
            mantissa_bits_within_error_bound = 0;
          }

          offset_bits = mantissa_bits_within_error_bound;

          if(offset_bits > 0) //has mantissa bits
          {
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            } 
          }
          else //no mantissa bit
          {
            decompressed_num++;
            decompressed[decompressed_num-1] = decompress_bitwise_double(bits, bits_num, before_value1, before_value2, before_value3);
            // printf("%f ", decompressed[decompressed_num-1]);

            if(before_value3 == -1) 
            {
              before_value3 = decompressed[decompressed_num-1]; 
            }
            else if(before_value2 == -1) 
            {
              before_value2 = decompressed[decompressed_num-1];
            }
            else if(before_value1 == -1) 
            {
              before_value1 = decompressed[decompressed_num-1];
            }
            else
            {
              before_value3 = before_value2;
              before_value2 = before_value1;
              before_value1 = decompressed[decompressed_num-1];
            }

            bits = NULL;
            bits_num = 0;

            if(bit == 0)
            {
              offset_bits = 1+11;
              bits_num++;
              bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
              if (bits_more != NULL) 
              {
                bits = bits_more;
                bits[bits_num-1] = bit + '0';
              }
              else 
              {
                free(bits);
                printf("Error (re)allocating memory\n");
                exit(1);
              }             
            }
            else if(bit == 1)
            {
              offset_bits = 3;
              bits_num++;
              bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
              if (bits_more != NULL) 
              {
                bits = bits_more;
                bits[bits_num-1] = bit + '0';
              }
              else 
              {
                free(bits);
                printf("Error (re)allocating memory\n");
                exit(1);
              }             
            }              
          }
        }
      }
      else
      {
        bits_num++;
        bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
        if (bits_more != NULL) 
        {
          bits = bits_more;
          bits[bits_num-1] = bit + '0';
        }
        else 
        {
          free(bits);
          printf("Error (re)allocating memory\n");
          exit(1);
        }        
      }
      offset_bits--;
      if(offset_bits == 0 && bits_num != 1+11)
      {
        decompressed_num++;
        decompressed[decompressed_num-1] = decompress_bitwise_double(bits, bits_num, before_value1, before_value2, before_value3);
        // printf("%f ", decompressed[decompressed_num-1]);
        
        if(before_value3 == -1) 
        {
          before_value3 = decompressed[decompressed_num-1]; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = decompressed[decompressed_num-1];
        }
        else if(before_value1 == -1) 
        {
          before_value1 = decompressed[decompressed_num-1];
        }
        else
        {
          before_value3 = before_value2;
          before_value2 = before_value1;
          before_value1 = decompressed[decompressed_num-1];
        }

        bits = NULL;
        bits_num = 0;
      }
    }       
  }
  return decompressed;
}

double decompress_bitwise_double(char* bits, int bits_num, double before_value1, double before_value2, double before_value3)
{
  if(bits_num == 3)
  {
    if(bits[0] == '1')
    {
      if(bits[1] == '0' && bits[2] == '0')
      {
        return 0.0;
      }
      else if(bits[1] == '0' && bits[2] == '1')
      {
        return before_value1;
      }
      else if(bits[1] == '1' && bits[2] == '0')
      {
        return 2*before_value1 - before_value2;
      }
      else if(bits[1] == '1' && bits[2] == '1')
      {
        return 3*before_value1 - 3*before_value2 + before_value3;
      }
    }
    else
    {
      printf("Error start bit of 3 bits is 0\n");
      exit(1);
    }
  }
  else
  {
    if(bits_num == sizeof(double)*8)
    {
      return strtodbl(bits);
    }
    else
    {
      char* bits64 = (char*)realloc(bits, sizeof(double)*8);
      bits64[bits_num] = '1';
      if(bits_num+1 < sizeof(double)*8)     
      {
        for(int i=bits_num+1; i< sizeof(double)*8; i++)
        {
          bits64[i] = '0';
        }
      }
      return strtodbl(bits64); 
    }
  }
}

float* myDecompress_bitwise(unsigned char* data_bits, int bytes, int num)
{
  int offset_bits = 0;
  char* bits = NULL;
  char* bits_more = NULL;
  int bits_num = 0;
  int min_shift = 0;

  float before_value1=-1, before_value2=-1, before_value3=-1;
  float* decompressed = (float*) malloc(sizeof(float)*num);
  int decompressed_num = 0;

  for(int i=0; i<bytes; i++)
  {
    //if(i == bytes - 1 && pos != 8) min_shift = pos;

    for (int j=7; j>=min_shift; j--) //each bit of byte
    {
      int bit = (data_bits[i] >> j) & 1;

      //printf("%d(%d)", bit, bits_num);

      if(offset_bits == 0) //start bit
      {
        if(bits_num == 0) //not start bit of mantissa
        {
          if(bit == 0)
          {
            offset_bits = 1+8;
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            }             
          }
          else if(bit == 1)
          {
            offset_bits = 3; //100, 101, 110, 111
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            }             
          }
        }
        else //start bit of mantissa
        {
          int expo_value = 0;
          for(int i=1; i<9; i++)
          {
            expo_value += (bits[i]-'0')*pow(2,8-i);
          }
          expo_value -= 127;           

          if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

          int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
          if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
          {
            mantissa_bits_within_error_bound = 23;
          }
          else if(mantissa_bits_within_error_bound < 0)
          {
            mantissa_bits_within_error_bound = 0;
          }

          offset_bits = mantissa_bits_within_error_bound;

          if(offset_bits > 0) //has mantissa bits
          {
            bits_num++;
            bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
            if (bits_more != NULL) 
            {
              bits = bits_more;
              bits[bits_num-1] = bit + '0';
            }
            else 
            {
              free(bits);
              printf("Error (re)allocating memory\n");
              exit(1);
            } 
          }
          else //no mantissa bit
          {
            decompressed_num++;
            decompressed[decompressed_num-1] = decompress_bitwise_float(bits, bits_num, before_value1, before_value2, before_value3);
            // printf("%f ", decompressed[decompressed_num-1]);

            if(before_value3 == -1) 
            {
              before_value3 = decompressed[decompressed_num-1]; 
            }
            else if(before_value2 == -1) 
            {
              before_value2 = decompressed[decompressed_num-1];
            }
            else if(before_value1 == -1) 
            {
              before_value1 = decompressed[decompressed_num-1];
            }
            else
            {
              before_value3 = before_value2;
              before_value2 = before_value1;
              before_value1 = decompressed[decompressed_num-1];
            }

            bits = NULL;
            bits_num = 0;

            if(bit == 0)
            {
              offset_bits = 1+8;
              bits_num++;
              bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
              if (bits_more != NULL) 
              {
                bits = bits_more;
                bits[bits_num-1] = bit + '0';
              }
              else 
              {
                free(bits);
                printf("Error (re)allocating memory\n");
                exit(1);
              }             
            }
            else if(bit == 1)
            {
              offset_bits = 3;
              bits_num++;
              bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
              if (bits_more != NULL) 
              {
                bits = bits_more;
                bits[bits_num-1] = bit + '0';
              }
              else 
              {
                free(bits);
                printf("Error (re)allocating memory\n");
                exit(1);
              }             
            }              
          }
        }
      }
      else
      {
        bits_num++;
        bits_more = (char*)realloc(bits, sizeof(char)*bits_num);
        if (bits_more != NULL) 
        {
          bits = bits_more;
          bits[bits_num-1] = bit + '0';
        }
        else 
        {
          free(bits);
          printf("Error (re)allocating memory\n");
          exit(1);
        }        
      }
      offset_bits--;
      if(offset_bits == 0 && bits_num != 1+8)
      {
        decompressed_num++;
        decompressed[decompressed_num-1] = decompress_bitwise_float(bits, bits_num, before_value1, before_value2, before_value3);
        // printf("%f ", decompressed[decompressed_num-1]);
        
        if(before_value3 == -1) 
        {
          before_value3 = decompressed[decompressed_num-1]; 
        }
        else if(before_value2 == -1) 
        {
          before_value2 = decompressed[decompressed_num-1];
        }
        else if(before_value1 == -1) 
        {
          before_value1 = decompressed[decompressed_num-1];
        }
        else
        {
          before_value3 = before_value2;
          before_value2 = before_value1;
          before_value1 = decompressed[decompressed_num-1];
        }

        bits = NULL;
        bits_num = 0;
      }
    }       
  }
  return decompressed;
}

float decompress_bitwise_float(char* bits, int bits_num, float before_value1, float before_value2, float before_value3)
{
  if(bits_num == 3)
  {
    if(bits[0] == '1')
    {
      if(bits[1] == '0' && bits[2] == '0')
      {
        return 0.0;
      }
      else if(bits[1] == '0' && bits[2] == '1')
      {
        return before_value1;
      }
      else if(bits[1] == '1' && bits[2] == '0')
      {
        return 2*before_value1 - before_value2;
      }
      else if(bits[1] == '1' && bits[2] == '1')
      {
        return 3*before_value1 - 3*before_value2 + before_value3;
      }
    }
    else
    {
      printf("Error start bit of 3 bits is 0\n");
      exit(1);
    }
  }
  else
  {
    if(bits_num == sizeof(float)*8)
    {
      return strtofloat(bits);
    }
    else
    {
      char* bits32 = (char*)realloc(bits, sizeof(float)*8);
      bits32[bits_num] = '1';
      if(bits_num+1 < sizeof(float)*8)     
      {
        for(int i=bits_num+1; i< sizeof(float)*8; i++)
        {
          bits32[i] = '0';
        }
      }
      return strtofloat(bits32); 
    }
  }
}

//bitwise myCompress for k-means (double)
void myCompress_bitwise_double(double data[], int num, unsigned char** data_bits, int* bytes, int* pos)
{
  double real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  double diff1, diff2, diff3, diff_min, selected_predict_value;
  char compress_type;
  int a=0, b=0, c=0, d=0;

  // unsigned char* data_bits = NULL;
  // int flag = 0; //0, 1
  // int bytes = 0; //total bytes of compressed data
  // int pos = 8; //position of filled bit in last byte --> 87654321

  for(int n=0; n<num; n++)
  {
    real_value = data[n];

    if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    {
      //if(real_value == 0)
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
      }
      else
      {
        compress_bitwise_double(real_value, data_bits, bytes, pos);
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
      compress_type = 'a'; //101
      selected_predict_value = predict_value1;
      if(diff2<diff_min)
      {
        diff_min = diff2;
        compress_type = 'b'; //110
        selected_predict_value = predict_value2;
      }
      if(diff3<diff_min)
      {
        diff_min = diff3;
        compress_type = 'c'; //111
        selected_predict_value = predict_value3;
      }        

      before_value3 = before_value2;
      before_value2 = before_value1;
      before_value1 = real_value;
      
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
      }
      else if(diff_min<=absErrorBound) 
      {
        if(compress_type == 'a')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);
          add_bit_to_bytes(data_bits, bytes, pos, 1);    
          a++;    
        }
        else if(compress_type == 'b')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);  
          b++;
        }
        else if(compress_type == 'c')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);  
          c++;
        }
        else
        {
          printf("Error compress_type\n");
          exit(1);
        }
      }
      else 
      {
        compress_bitwise_double(real_value, data_bits, bytes, pos);            
      }
    }
  }

  //printf("compression pattern: a = %d (%f), b = %d (%f), c = %d (%f), d = %d (%f), num = %d\n", a, (float)a/num, b, (float)b/num, c, (float)c/num, d, (float)d/num, num);
}

//bitwise myCompress for ping-pong & himeno (float)
void myCompress_bitwise(float data[], int num, unsigned char** data_bits, int* bytes, int* pos)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  float diff1, diff2, diff3, diff_min, selected_predict_value;
  char compress_type;
  int a=0, b=0, c=0, d=0;

  // unsigned char* data_bits = NULL;
  // int flag = 0; //0, 1
  // int bytes = 0; //total bytes of compressed data
  // int pos = 8; //position of filled bit in last byte --> 87654321

  // FILE *fp;
  // fp = fopen("bitcomp.txt", "w");
  // assert(fp);

  for(int n=0; n<num; n++)
  {
    real_value = data[n];

    if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    {
      //if(real_value == 0)
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
        //printf("100\n");
        // fprintf(fp, "1 0 0\n");
      }
      else
      {
        compress_bitwise_float(real_value, data_bits, bytes, pos);
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
      compress_type = 'a'; //101
      selected_predict_value = predict_value1;
      if(diff2<diff_min)
      {
        diff_min = diff2;
        compress_type = 'b'; //110
        selected_predict_value = predict_value2;
      }
      if(diff3<diff_min)
      {
        diff_min = diff3;
        compress_type = 'c'; //111
        selected_predict_value = predict_value3;
      }        

      before_value3 = before_value2;
      before_value2 = before_value1;
      before_value1 = real_value;
      
      if(fabs(real_value) < absErrorBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
        //printf("100\n");
        // fprintf(fp, "1 0 0\n");
      }
      else if(diff_min<=absErrorBound) 
      {
        if(compress_type == 'a')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);
          add_bit_to_bytes(data_bits, bytes, pos, 1);   
          a++;     
          //printf("101\n");
          // fprintf(fp, "1 0 1\n");
        }
        else if(compress_type == 'b')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 0);  
          b++;
          //printf("110\n");
          // fprintf(fp, "1 1 0\n");
        }
        else if(compress_type == 'c')
        {
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);
          add_bit_to_bytes(data_bits, bytes, pos, 1);  
          c++;
          //printf("111\n");
          // fprintf(fp, "1 1 1\n");
        }
        else
        {
          printf("Error compress_type\n");
          exit(1);
        }
      }
      else 
      {
        compress_bitwise_float(real_value, data_bits, bytes, pos);            
      }
    }
  }

  // fclose(fp);

  //printf("compression pattern: a = %d (%f), b = %d (%f), c = %d (%f), d = %d (%f), num = %d\n", a, (float)a/num, b, (float)b/num, c, (float)c/num, d, (float)d/num, num);
}

void compress_bitwise_double(double real_value, unsigned char** data_bits, int* bytes, int* pos)
{
  double double10 = real_value;
  char double_arr[64+1];
  doubletostr(&double10, double_arr);

  int expo_value = 0;
  for(int i=1; i<12; i++)
  {
    expo_value += (double_arr[i]-'0')*pow(2,11-i);
  }
  expo_value -= 1023;  

  if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound); 

  int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;

  if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
  {
    mantissa_bits_within_error_bound = 52;
  }
  else if(mantissa_bits_within_error_bound < 0)
  {
    mantissa_bits_within_error_bound = 0;
  }
  int bits_after_compress = 1+11+mantissa_bits_within_error_bound;  

  for(int i=0; i<bits_after_compress; i++)
  {
    add_bit_to_bytes(data_bits, bytes, pos, double_arr[i]-'0');
  }
}

void compress_bitwise_float(float real_value, unsigned char** data_bits, int* bytes, int* pos)
{
  float float10 = real_value;
  char float_arr[32+1];
  floattostr(&float10, float_arr);

  int expo_value = 0;
  for(int i=1; i<9; i++)
  {
    expo_value += (float_arr[i]-'0')*pow(2,8-i);
  }
  expo_value -= 127;  

  if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

  int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;

  if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
  {
    mantissa_bits_within_error_bound = 23;
  }
  else if(mantissa_bits_within_error_bound < 0)
  {
    mantissa_bits_within_error_bound = 0;
  }

  int bits_after_compress = 1+8+mantissa_bits_within_error_bound;  

  // FILE *fp;
  // fp = fopen("bitcomp.txt", "a");
  // assert(fp);

  for(int i=0; i<bits_after_compress; i++)
  {
    add_bit_to_bytes(data_bits, bytes, pos, float_arr[i]-'0');
    //printf("%d", float_arr[i]-'0');
    // fprintf(fp, "%d ", float_arr[i]-'0');
  }
  //printf("\n");
  // fprintf(fp, "\n");
  // fclose(fp);
}

double toSmallDataset_double(double data[], double** data_small, int num)
{
  *data_small = malloc(sizeof(double) * num);
  double min = data[0];

  for(int i=1; i<num; i++)
  {
    if(data[i]<min)
    {
      min = data[i];
    }
  }

  for(int i=0; i<num; i++)
  {
    (*data_small)[i] = data[i] - min;
  }

  return min;
}

float toSmallDataset_float(float data[], float** data_small, int num)
{
  *data_small = malloc(sizeof(float) * num);
  float min = data[0];

  for(int i=1; i<num; i++)
  {
    if(data[i]<min)
    {
      min = data[i];
    }
  }

  for(int i=0; i<num; i++)
  {
    (*data_small)[i] = data[i] - min;
  }

  return min;
}

double med_dataset_double(double* data, int num, int* type)
{
  double total = 0;
  double max = data[0];
  for(int i=0; i<num; i++)
  {
    total += data[i];
    if(data[i] > max)
    {
      max = data[i];
    }
  }
  int add = 0;
  for(int i=10; i>0; i--)
  {
    add += pow(2, i);
    if(max < add - 1023)
    {
      *type = 11 - i;
      break;
    }
  }  
  return total/num;
  // double medium = total/num; //min = 0
  // if(medium > max/2) return medium + (max - medium)/2;
  // else if(medium < max/2) return medium/2;
  // else return medium;
}

float med_dataset_float(float* data, int num, int* type)
{
  float total = 0;
  float max = data[0];
  for(int i=0; i<num; i++)
  {
    total += data[i];
    if(data[i] > max)
    {
      max = data[i];
    }
  }
  int add = 0;
  for(int i=7; i>0; i--)
  {
    add += pow(2, i);
    if(max < add - 127)
    {
      *type = 8 - i;
      break;
    }
  }
  return total/num;
  // float medium = total/num; //min = 0
  // if(medium > max/2) return medium + (max - medium)/2;
  // else if(medium < max/2) return medium/2;
  // else return medium;
}

float calCompressRatio_bitwise_double2(float data[], int num)
{
  int bits_after_compress = 0;

  for(int n=0; n<num; n++)
  {
    double double10 = data[n];
    char double_arr[64+1];
    doubletostr(&double10, double_arr);
    //printf("%s \n", double_arr); 

    int expo_value = 0;
    //printf("%f \n", double10); 
    for(int i=1; i<12; i++)
    {
      expo_value += (double_arr[i]-'0')*pow(2,11-i);
      //printf("%d ", double_arr[i]-'0'); 
    }
    expo_value -= 1023;

    //printf("%d ", expo_value);        

    if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

    int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;

    if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
    {
      mantissa_bits_within_error_bound = 52;
    }
    else if(mantissa_bits_within_error_bound < 0)
    {
      mantissa_bits_within_error_bound = 0;
    }
    bits_after_compress += 1+11+mantissa_bits_within_error_bound;  
  }

  return (float)bits_after_compress/(sizeof(double)*8*num);
}

float calCompressRatio_bitwise_double(double data[], int num)
{
  int bits_after_compress = 0;

  for(int n=0; n<num; n++)
  {
    double double10 = data[n];
    char double_arr[64+1];
    doubletostr(&double10, double_arr);
    //printf("%s \n", double_arr); 

    int expo_value = 0;
    //printf("%f \n", double10); 
    for(int i=1; i<12; i++)
    {
      expo_value += (double_arr[i]-'0')*pow(2,11-i);
      //printf("%d ", double_arr[i]-'0'); 
    }
    expo_value -= 1023;

    //printf("%d ", expo_value);  

    if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);     

    int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;

    if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
    {
      mantissa_bits_within_error_bound = 52;
    }
    else if(mantissa_bits_within_error_bound < 0)
    {
      mantissa_bits_within_error_bound = 0;
    }
    bits_after_compress += 1+11+mantissa_bits_within_error_bound;  
  }

  return (float)bits_after_compress/(sizeof(double)*8*num);
}

float calCompressRatio_bitwise_float(float data[], int num)
{
  int bits_after_compress = 0;

  for(int n=0; n<num; n++)
  {
    float float10 = data[n];
    char float_arr[32+1];
    floattostr(&float10, float_arr);
    //printf("%s \n", float_arr); 

    int expo_value = 0;
    for(int i=1; i<9; i++)
    {
      expo_value += (float_arr[i]-'0')*pow(2,8-i);
      //printf("%d ", float_arr[i]-'0'); 
    }
    expo_value -= 127;

    //printf("%d ", expo_value);   

    if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

    int mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;

    if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
    {
      mantissa_bits_within_error_bound = 23;
    }
    else if(mantissa_bits_within_error_bound < 0)
    {
      mantissa_bits_within_error_bound = 0;
    }
    bits_after_compress += 1+8+mantissa_bits_within_error_bound;  
  }

  return (float)bits_after_compress/(sizeof(float)*8*num);
}

float* transform_3d_array_to_1d_array(float data[MIMAX][MJMAX][MKMAX], int ijk, int v, int imax, int jmax, int kmax)
{
  int A, B;
  //float array_1d[len]; 
  float* array_1d;
  
  if(ijk == 1) 
  {
    A = jmax;
    B = kmax;
  }
  else if(ijk == 2) 
  {
    A = imax;
    B = kmax;
  }
  else if(ijk == 3)
  {
    A = imax;
    B = jmax;
  } 
  array_1d = (float*) malloc(sizeof(float)*A*B);

  int array_1d_len = 0;
  for(int a=0; a<A; a++)
  {
    for(int b=0; b<B; b++)
    {
      if(ijk == 1) array_1d[array_1d_len++] = data[v][a][b]; // v is const
      else if(ijk == 2) array_1d[array_1d_len++] = data[a][v][b];
      else if(ijk == 3) array_1d[array_1d_len++] = data[a][b][v];
    }
  } 
  return array_1d;  
}

//myDecompress for k-means (double)
double* myDecompress_double(double array_double[], char array_char[], int array_char_displacement[], int num)
{
  double* data = (double*) malloc(sizeof(double)*num);
  int array_double_p = 0, array_char_p = 0, array_char_displacement_p = 0;
  for(int i=0; i<num; i++)
  {
    if(array_char_displacement[array_char_displacement_p] - 1 == i)
    {
      if(array_char[array_char_p] == 'a')
      {
        data[i] = data[i-1];
      }
      else if(array_char[array_char_p] == 'b')
      {
        data[i] = 2*data[i-1] - data[i-2];
      }
      else if(array_char[array_char_p] == 'c')
      {
        data[i] = 3*data[i-1] - 3*data[i-2] + data[i-3];
      }
      else if(array_char[array_char_p] == 'd')
      {
        data[i] = 4*data[i-1] - 6*data[i-2] + 4*data[i-3] - data[i-4];
      }      
      array_char_p++;
      array_char_displacement_p++;
    }
    else
    {
      data[i] = array_double[array_double_p];
      array_double_p++;
    }
  }
  return data;
}

//myCompress for k-means (double)
int myCompress_double(double data[], double** array_double, char** array_char, int** array_char_displacement, int num)
{
  double real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value1, predict_value2, predict_value3, predict_value4;
  double diff1, diff2, diff3, diff4, diff_min, selected_predict_value;
  int array_double_len = 0, array_char_len = 0;
  char compress_type;
  double* array_double_more = NULL;
  char* array_char_more = NULL;
  int* array_char_displacement_more = NULL;

  for(int n=0; n<num; n++)
  {
    real_value = data[n];

    if(before_value4 == -1 || before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    {
      array_double_len++;
      array_double_more = (double*)realloc(*array_double, sizeof(double)*array_double_len);
      if (array_double_more != NULL) 
      {
        *array_double = array_double_more;
        (*array_double)[array_double_len-1] = real_value;
      }
      else 
      {
        free(*array_double);
        printf("Error (re)allocating memory\n");
        exit(1);
      }        

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
      predict_value1 = before_value1;
      predict_value2 = 2*before_value1 - before_value2;
      predict_value3 = 3*before_value1 - 3*before_value2 + before_value3;
      predict_value4 = 4*before_value1 - 6*before_value2 + 4*before_value3 - before_value4;

      diff1 = fabs(predict_value1-real_value);
      diff2 = fabs(predict_value2-real_value);
      diff3 = fabs(predict_value3-real_value);
      diff4 = fabs(predict_value4-real_value);

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
      if(diff4<diff_min)
      {
        diff_min = diff4;
        compress_type = 'd';
        selected_predict_value = predict_value4;
      }             

      before_value4 = before_value3;
      before_value3 = before_value2;
      before_value2 = before_value1;
      before_value1 = real_value;
      
      if(diff_min<=absErrorBound) 
      {
        array_char_len++;
        array_char_more = (char*)realloc(*array_char, sizeof(char)*array_char_len);
        array_char_displacement_more = (int*)realloc(*array_char_displacement, sizeof(int)*array_char_len);
        if (array_char_more != NULL && array_char_displacement_more != NULL) 
        {
          *array_char = array_char_more;
          (*array_char)[array_char_len-1] = compress_type;
          *array_char_displacement = array_char_displacement_more;
          (*array_char_displacement)[array_char_len-1] = array_double_len + array_char_len;
        }
        else 
        {
          free(*array_char);
          free(*array_char_displacement);
          printf("Error (re)allocating memory\n");
          exit(1);
        } 
      }
      else 
      {
        array_double_len++;
        array_double_more = (double*)realloc(*array_double, sizeof(double)*array_double_len);
        if (array_double_more != NULL) 
        {
          *array_double = array_double_more;
          (*array_double)[array_double_len-1] = real_value;
        }
        else 
        {
          free(*array_double);
          printf("Error (re)allocating memory\n");
          exit(1);
        }             
      }
    }
  }
  return array_double_len;
}

//myDecompress for ping-pong & himeno (float)
float* myDecompress(float array_float[], char array_char[], int array_char_displacement[], int num)
{
  float* data = (float*) malloc(sizeof(float)*num);
  int array_float_p = 0, array_char_p = 0, array_char_displacement_p = 0;
  for(int i=0; i<num; i++)
  {
    if(array_char_displacement[array_char_displacement_p] - 1 == i)
    {
      if(array_char[array_char_p] == 'a')
      {
        data[i] = data[i-1];
      }
      else if(array_char[array_char_p] == 'b')
      {
        data[i] = 2*data[i-1] - data[i-2];
      }
      else if(array_char[array_char_p] == 'c')
      {
        data[i] = 3*data[i-1] - 3*data[i-2] + data[i-3];
      }
      else if(array_char[array_char_p] == 'd')
      {
        data[i] = 4*data[i-1] - 6*data[i-2] + 4*data[i-3] - data[i-4];
      }      
      array_char_p++;
      array_char_displacement_p++;
    }
    else
    {
      data[i] = array_float[array_float_p];
      array_float_p++;
    }
  }
  return data;
}

//myCompress for ping-pong & himeno (float)
int myCompress(float data[], float** array_float, char** array_char, int** array_char_displacement, int num)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value1, predict_value2, predict_value3, predict_value4;
  float diff1, diff2, diff3, diff4, diff_min, selected_predict_value;
  int array_float_len = 0, array_char_len = 0;
  char compress_type;
  //float compress_ratio;
  // float* array_float = NULL; //(float*)malloc(sizeof(float));
  float* array_float_more = NULL;
  // char* array_char = NULL; //(char*)malloc(sizeof(char));
  char* array_char_more = NULL;
  // int* array_char_displacement = NULL;
  int* array_char_displacement_more = NULL;

  for(int n=0; n<num; n++)
  {
    real_value = data[n];

    if(before_value4 == -1 || before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    {
      array_float_len++;
      array_float_more = (float*)realloc(*array_float, sizeof(float)*array_float_len);
      if (array_float_more != NULL) 
      {
        *array_float = array_float_more;
        (*array_float)[array_float_len-1] = real_value;
      }
      else 
      {
        free(*array_float);
        printf("Error (re)allocating memory\n");
        exit(1);
      }        

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
      predict_value1 = before_value1;
      predict_value2 = 2*before_value1 - before_value2;
      predict_value3 = 3*before_value1 - 3*before_value2 + before_value3;
      predict_value4 = 4*before_value1 - 6*before_value2 + 4*before_value3 - before_value4;

      diff1 = fabs(predict_value1-real_value);
      diff2 = fabs(predict_value2-real_value);
      diff3 = fabs(predict_value3-real_value);
      diff4 = fabs(predict_value4-real_value);

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
      if(diff4<diff_min)
      {
        diff_min = diff4;
        compress_type = 'd';
        selected_predict_value = predict_value4;
      } 

      before_value4 = before_value3;
      before_value3 = before_value2;
      before_value2 = before_value1;
      before_value1 = real_value;
      
      if(diff_min<=absErrorBound) 
      {
        array_char_len++;
        array_char_more = (char*)realloc(*array_char, sizeof(char)*array_char_len);
        array_char_displacement_more = (int*)realloc(*array_char_displacement, sizeof(int)*array_char_len);
        if (array_char_more != NULL && array_char_displacement_more != NULL) 
        {
          *array_char = array_char_more;
          (*array_char)[array_char_len-1] = compress_type;
          *array_char_displacement = array_char_displacement_more;
          (*array_char_displacement)[array_char_len-1] = array_float_len + array_char_len;
        }
        else 
        {
          free(*array_char);
          free(*array_char_displacement);
          printf("Error (re)allocating memory\n");
          exit(1);
        } 
      }
      else 
      {
        array_float_len++;
        array_float_more = (float*)realloc(*array_float, sizeof(float)*array_float_len);
        if (array_float_more != NULL) 
        {
          *array_float = array_float_more;
          (*array_float)[array_float_len-1] = real_value;
        }
        else 
        {
          free(*array_float);
          printf("Error (re)allocating memory\n");
          exit(1);
        }             
      }
    }
  }
  // if(byte_or_bit == 1)
  // {
  //   compress_ratio = (float)(array_char_len*sizeof(char)+array_float_len*sizeof(float))/((array_char_len+array_float_len)*sizeof(float));
  // }
  // else if(byte_or_bit == 2)
  // {
  //   compress_ratio = (float)(array_char_len*2+array_float_len*sizeof(float)*8)/((array_char_len+array_float_len)*sizeof(float)*8);
  // }  
  // printf("Compression rate: %f \n", 1/compress_ratio);
  return array_float_len;
}

//myCompress for himeno
float calcCompressionRatio_himeno_ij_ik_jk(float data[MIMAX][MJMAX][MKMAX], int ijk, int v, int imax, int jmax, int kmax)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value1, predict_value2, predict_value3, predict_value4;
  float diff1, diff2, diff3, diff4, diff_min, selected_predict_value;
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
    A = jmax;
    B = kmax;
  }
  else if(ijk == 2) 
  {
    A = imax;
    B = kmax;
  }
  else if(ijk == 3)
  {
    A = imax;
    B = jmax;
  } 
  for(int a=0; a<A; a++)
  {
    for(int b=0; b<B; b++)
    {
      if(ijk == 1) real_value = data[v][a][b]; // v is const
      else if(ijk == 2) real_value = data[a][v][b];
      else if(ijk == 3) real_value = data[a][b][v];

      if(before_value4 == -1 || before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
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
          printf("Error (re)allocating memory\n");
          exit(1);
        }        

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
        predict_value1 = before_value1;
        predict_value2 = 2*before_value1 - before_value2;
        predict_value3 = 3*before_value1 - 3*before_value2 + before_value3;
        predict_value4 = 4*before_value1 - 6*before_value2 + 4*before_value3 - before_value4;

        diff1 = fabs(predict_value1-real_value);
        diff2 = fabs(predict_value2-real_value);
        diff3 = fabs(predict_value3-real_value);
        diff4 = fabs(predict_value4-real_value);

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
        if(diff4<diff_min)
        {
          diff_min = diff4;
          compress_type = 'd';
          selected_predict_value = predict_value4;
        }              

        before_value4 = before_value3;
        before_value3 = before_value2;
        before_value2 = before_value1;
        before_value1 = real_value;

        if(diff_min<=absErrorBound) 
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
            printf("Error (re)allocating memory\n");
            exit(1);
          } 
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
            printf("Error (re)allocating memory\n");
            exit(1);
          }             
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

float calcCompressionRatio_himeno_sz(float data[MIMAX][MJMAX][MKMAX], int ijk, int v, int imax, int jmax, int kmax)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  float diff1, diff2, diff3, diff_min, predict_diff, selected_predict_value;
  char compress_type;
  float compress_ratio;
  int A, B;
  long origin_bits=0, compressed_bits=0;

  if(ijk == 1) 
  {
    A = jmax;
    B = kmax;
  }
  else if(ijk == 2) 
  {
    A = imax;
    B = kmax;
  }
  else if(ijk == 3)
  {
    A = imax;
    B = jmax;
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
        before_value1 = real_value;

        if(diff_min<=absErrorBound) 
        {
          if(byte_or_bit == 1)
          {
            compressed_bits += sizeof(char)*8; 
          }
          else if(byte_or_bit == 2)
          {
            compressed_bits += 2; 
          }
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

          if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

          mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
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
            if(mantissa_bits_within_error_bound%8 != 0) compressed_bits += 1+8+(mantissa_bits_within_error_bound/8+1)*8;  
            else compressed_bits += 1+8+mantissa_bits_within_error_bound;
          }
          else if(byte_or_bit == 2)
          {
            compressed_bits += 1+8+mantissa_bits_within_error_bound;  
          }
        }
      }
    }
  } 
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_himeno_nolossy_performance(float data[MIMAX][MJMAX][MKMAX], int ijk, int v, int imax, int jmax, int kmax)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  float diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;
  int A, B;

  if(ijk == 1) 
  {
    A = jmax;
    B = kmax;
  }
  else if(ijk == 2) 
  {
    A = imax;
    B = kmax;
  }
  else if(ijk == 3)
  {
    A = imax;
    B = jmax;
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
        before_value1 = real_value;

        char c[sizeof(float)*8];
        getFloatBin(diff4, c);
        for(int i=1;i<sizeof(float)*8;i++)
        {
          if(c[i] != 0) 
          {
            if(byte_or_bit == 1)
            {
              if((sizeof(float)*8 - i + 3 + 1)%8 != 0) compressed_bits += ((sizeof(float)*8 - i + 3 + 1)/8+1)*8;  
              else compressed_bits += sizeof(float)*8 - i + 3 + 1;
            }
            else if(byte_or_bit == 2)
            {
              compressed_bits += sizeof(float)*8 - i + 3 + 1;
            }
            break;
          } 
        }
      }
    }
  }
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_himeno_nolossy_area(float data[MIMAX][MJMAX][MKMAX], int ijk, int v, int imax, int jmax, int kmax)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  float diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;
  int A, B;
  int cdb = 512, cdb_num = 1, occupied_bits = 0, indication = 5, le = 285, re = 222, re1 = 2, re2 = 4, re3 = 32, llrb = 2, ex = 1;

  if(ijk == 1) 
  {
    A = jmax;
    B = kmax;
  }
  else if(ijk == 2) 
  {
    A = imax;
    B = kmax;
  }
  else if(ijk == 3)
  {
    A = imax;
    B = jmax;
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
        before_value1 = real_value;

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
      }
    }
  }
  compressed_bits = cdb_num*cdb;
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_sz_float(float data[], int num)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  float diff1, diff2, diff3, diff_min, predict_diff, selected_predict_value;
  char compress_type;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;

  for(int n=0; n<num; n++)
  {
    real_value = data[n]; 

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
      before_value1 = real_value;

      if(diff_min<=absErrorBound) 
      {
        if(byte_or_bit == 1)
        {
          compressed_bits += sizeof(char)*8; 
        }
        else if(byte_or_bit == 2)
        {
          compressed_bits += 2; 
        }
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

        if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

        mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
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
          if(mantissa_bits_within_error_bound%8 != 0) compressed_bits += 1+8+(mantissa_bits_within_error_bound/8+1)*8;  
          else compressed_bits += 1+8+mantissa_bits_within_error_bound;
        }
        else if(byte_or_bit == 2)
        {
          compressed_bits += 1+8+mantissa_bits_within_error_bound;  
        }
      }
    }
  }
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_nolossy_performance_float(float data[], int num)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  float diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;

  for(int n=0; n<num; n++)
  {
    real_value = data[n]; 

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
      before_value1 = real_value;

      char c[sizeof(float)*8];
      getFloatBin(diff4, c);
      for(int i=1;i<sizeof(float)*8;i++)
      {
        if(c[i] != 0) 
        {
          if(byte_or_bit == 1)
          {
            if((sizeof(float)*8 - i + 3 + 1)%8 != 0) compressed_bits += ((sizeof(float)*8 - i + 3 + 1)/8+1)*8;  
            else compressed_bits += sizeof(float)*8 - i + 3 + 1;
          }
          else if(byte_or_bit == 2)
          {
            compressed_bits += sizeof(float)*8 - i + 3 + 1;
          }
          break;
        } 
      }
    }    
  }
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_nolossy_area_float(float data[], int num)
{
  float real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  float diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;
  int cdb = 512, cdb_num = 1, occupied_bits = 0, indication = 5, le = 285, re = 222, re1 = 2, re2 = 4, re3 = 32, llrb = 2, ex = 1;

  for(int n=0; n<num; n++)
  {
    real_value = data[n]; 

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
      before_value1 = real_value;

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
    }    
  }
  compressed_bits = cdb_num*cdb;
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;  
}

float calcCompressionRatio_sz_double(double data[], int num)
{
  double real_value, before_value1=-1, before_value2=-1, before_value3=-1, predict_value1, predict_value2, predict_value3;
  double diff1, diff2, diff3, diff_min, predict_diff, selected_predict_value;
  char compress_type;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;

  for(int n=0; n<num; n++)
  {
    real_value = data[n]; 

    origin_bits += sizeof(double)*8;     

    if(before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    { 
      compressed_bits += sizeof(double)*8; 
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
      before_value1 = real_value;

      if(diff_min<=absErrorBound) 
      {
        if(byte_or_bit == 1)
        {
          compressed_bits += sizeof(char)*8; 
        }
        else if(byte_or_bit == 2)
        {
          compressed_bits += 2; 
        }
      }
      else 
      {
        double max, min;
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

        char c[sizeof(double)*8];
        getDoubleBin(predict_diff/2, c);
        int expo_value = 0;
        int mantissa_bits_within_error_bound;

        for(int i=1;i<12;i++) //1-9 exponential part of float (1-12 in the case of double)
        {
          if(c[i] != 0) 
          {
            expo_value += pow(2, 11-i);
          }  
        }
        expo_value -= 1023; 

        if(absErrorBound_binary == -100) absErrorBound_binary = to_absErrorBound_binary(absErrBound);

        mantissa_bits_within_error_bound = absErrorBound_binary + expo_value;
        if(mantissa_bits_within_error_bound > 52) //23 mantissa part of float (52 in the case of double)
        {
          mantissa_bits_within_error_bound = 52;
        }
        else if(mantissa_bits_within_error_bound < 0)
        {
          mantissa_bits_within_error_bound = 0;
        }
        if(byte_or_bit == 1)
        {
          if(mantissa_bits_within_error_bound%8 != 0) compressed_bits += 1+11+(mantissa_bits_within_error_bound/8+1)*8;  
          else compressed_bits += 1+11+mantissa_bits_within_error_bound;
        }
        else if(byte_or_bit == 2)
        {
          compressed_bits += 1+11+mantissa_bits_within_error_bound;  
        }
      }
    }
  }
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;
}

float calcCompressionRatio_nolossy_performance_double(double data[], int num)
{
  double real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  double diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;

  for(int n=0; n<num; n++)
  {
    real_value = data[n]; 

    origin_bits += sizeof(double)*8;

    if(before_value4 == -1 || before_value3 == -1 || before_value2 == -1 || before_value1 == -1)
    {
      compressed_bits += sizeof(double)*8;       
      
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
      before_value1 = real_value;

      char c[sizeof(double)*8];
      getDoubleBin(diff4, c);
      for(int i=1;i<sizeof(double)*8;i++)
      {
        if(c[i] != 0) 
        {
          if(byte_or_bit == 1)
          {
            if((sizeof(double)*8 - i + 3 + 1)%8 != 0) compressed_bits += ((sizeof(double)*8 - i + 3 + 1)/8+1)*8;  
            else compressed_bits += sizeof(double)*8 - i + 3 + 1;
          }
          else if(byte_or_bit == 2)
          {
            compressed_bits += sizeof(double)*8 - i + 3 + 1;
          }
          break;
        } 
      }
    }    
  }
  compress_ratio = (float)compressed_bits/origin_bits;
  return compress_ratio;  
}

float calcCompressionRatio_nolossy_area_double(double data[], int num)
{
  double real_value, before_value1=-1, before_value2=-1, before_value3=-1, before_value4=-1, predict_value4;
  double diff4;
  float compress_ratio;
  long origin_bits=0, compressed_bits=0;
  int cdb = 512, cdb_num = 1, occupied_bits = 0, indication = 5, le = 285, re = 222, re1 = 2, re2 = 4, re3 = 32, llrb = 2, ex = 1;

  for(int n=0; n<num; n++)
  {
    real_value = data[n]; 

    origin_bits += sizeof(double)*8;

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
      before_value1 = real_value;

      char c[sizeof(double)*8];
      getDoubleBin(diff4, c);
      for(int i=1;i<sizeof(double)*8;i++)
      {
        if(c[i] != 0) 
        {
          int nonzero = sizeof(double)*8 - i;
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

void getDoubleBin(double num,char bin[])
{
    int t = 1;
    int *f = (int*)(&num);
    for(int i=0;i<64;i++)
    {
        bin[i] = (*f)&(t<<63-i)?1:0;
    }
}

//100.0 --> 01000010110010000000000000000000
//str should have at least 33 byte.
void floattostr(float* a, char* str){
	unsigned int c;
	c= ((unsigned int*)a)[0]; 
	for(int i=0;i<32;i++){
		str[31-i]=(char)(c&1)+'0';
		c>>=1;
	}
	str[32] = '\0';
}

//100.0 --> 0100000001011001000000000000000000000000000000000000000000000000
//str should have at least 65 byte.
void doubletostr(double* a, char* str){
	long long c;
	c= ((long long*)a)[0]; 
	for(int i=0;i<64;i++){
		str[63-i]=(char)(c&1)+'0';
		c>>=1;
	}
	str[64] = '\0';
}

//01000010110010000000000000000000 --> 100.0
float strtofloat(char * str){
	unsigned int flt = 0;
	for(int i=0;i<31;i++){
		flt += (str[i]-'0');
		flt <<= 1;
	}
	flt += (str[31]-'0');
	float * ret = (float*)&flt;
	return *ret;
}

//0100000001011001000000000000000000000000000000000000000000000000 --> 100.0
double strtodbl(char * str){
	long long dbl = 0;
	for(int i=0;i<63;i++){
		dbl += (str[i]-'0');
		dbl <<= 1;
	}
	dbl +=(str[63]-'0');
	double* db = (double*)&dbl;
	return *db;
}     

void writetobinary_float(const char *file, float* data, int count)
{
    FILE *fp;
    fp = fopen(file, "wb");
    //fopen_s(&fp, file, "wb");
 
    if (fp == NULL)
    {
        printf("%sのオープンに失敗しました。\n", file);
        return;
    }

    fwrite(data, sizeof(float), count, fp);
    fclose(fp);
    printf("%sに保存しました。\n", file);
}

void writetobinary_double(const char *file, double* data, int count)
{
    FILE *fp;
    fp = fopen(file, "wb");
    //fopen_s(&fp, file, "wb");
 
    if (fp == NULL)
    {
        printf("%sのオープンに失敗しました。\n", file);
        return;
    }
    
    fwrite(data, sizeof(double), count, fp);
    fclose(fp);
    printf("%sに保存しました。\n", file);
}

void writetobinary_char(const char *file, unsigned char* data, int count)
{
    FILE *fp;
    fp = fopen(file, "wb");
    //fopen_s(&fp, file, "wb");
 
    if (fp == NULL)
    {
        printf("%sのオープンに失敗しました。\n", file);
        return;
    }
    
    fwrite(data, sizeof(char), count, fp);
    fclose(fp);
    printf("%sに保存しました。\n", file);
}

float* readfrombinary_float(const char *file, int count)
{
    FILE *fp;
    fp = fopen(file, "rb");
    //fopen_s(&fp, file, "rb");
 
    if (fp == NULL)
    {
        printf("%sのオープンに失敗しました。\n", file);
        exit(0);
    }
 
    //float arr[count];
    float *arr = (float*)malloc(sizeof(float) * count);
 
    fread(arr, sizeof(float), count, fp);
    fclose(fp);

    return arr;
}

double* readfrombinary_double(const char *file, int count)
{
    FILE *fp;
    fp = fopen(file, "rb");
    //fopen_s(&fp, file, "rb");
 
    if (fp == NULL)
    {
        printf("%sのオープンに失敗しました。\n", file);
        exit(0);
    }
 
    //double arr[count];
    double *arr = (double*)malloc(sizeof(double) * count);
 
    fread(arr, sizeof(double), count, fp);
    fclose(fp);

    return arr;
}

unsigned char* readfrombinary_char(const char *file, int* bytes_sz)
{
    FILE *fp;
    fp = fopen(file, "rb");
    if (fp == NULL)
    {
        printf("%sのオープンに失敗しました。\n", file);
        exit(0);
    }
 
    fseek(fp, 0, SEEK_END);
    *bytes_sz = ftell(fp);
    fclose(fp);

    fp = fopen(file, "rb");
    if (fp == NULL)
    {
        printf("%sのオープンに失敗しました。\n", file);
        exit(0);
    }

    unsigned char *arr = (unsigned char*)malloc(sizeof(char) * (*bytes_sz));
 
    fread(arr, sizeof(char), *bytes_sz, fp);
    fclose(fp);

    return arr;
}

float* readfrombinary_writetotxt_float(const char *binaryfile, const char *txtfile, int count)
{
    FILE *fp;
    fp = fopen(binaryfile, "rb");
    assert(fp);

    float *arr = malloc(sizeof(float) * count);
    fread(arr, sizeof(float), count, fp);
    fclose(fp);

    fp = fopen(txtfile, "w");
    assert(fp);
     
    for(int i=0; i<count; i++)
    {
      fprintf(fp, "%f\n", arr[i]);
    }
    fclose(fp);

    return arr;
}

double* readfrombinary_writetotxt_double(const char *binaryfile, const char *txtfile, int count)
{
    FILE *fp;
    fp = fopen(binaryfile, "rb");
    assert(fp);

    double *arr = malloc(sizeof(double) * count);
    fread(arr, sizeof(double), count, fp);
    fclose(fp);

    fp = fopen(txtfile, "w");
    assert(fp);

    for(int i=0; i<count; i++)
    {
      fprintf(fp, "%lf\n", arr[i]);
    }
    fclose(fp);

    return arr;
}

void add_bit_to_bytes(unsigned char** data_bits, int* bytes, int* pos, int flag)
{
  if(*pos > 0 && *pos < 9)
  {
    if(*pos == 8) 
    {
      (*bytes)++;
      unsigned char* data_bits_more = (unsigned char*)realloc(*data_bits, sizeof(char)*(*bytes));
      if (data_bits_more != NULL) 
      {
        *data_bits = data_bits_more;
        (*data_bits)[*bytes-1] = 0; //put all 8 bits to 0
        bit_set(&((*data_bits)[*bytes-1]), *pos, flag);
        (*pos)--;
      }
      else 
      {
        free(*data_bits);
        printf("Error (re)allocating memory\n");
        exit(1);
      }         
    }
    else{
      bit_set(&((*data_bits)[(*bytes)-1]), *pos, flag);
      (*pos)--;     
    }
    if(*pos == 0) *pos = 8;
  }
  else
  {
    printf("Error position value\n");
    return;
  }
}

// n*8 bits, position --> 87654321, flag --> 1, 0
void bit_set(unsigned char *p_data, unsigned char position, int flag)
{
	// int i = 0;
	assert(p_data);
	if (position > 8 || position < 1 || (flag != 0 && flag != 1))
	{
		printf("输入有误！\n");
		return;
	}
	if (flag != (*p_data >> (position - 1) & 1))
	{
		*p_data ^= 1 << (position - 1);
	}
	// for (i = 7; i >= 0; i--)     //由低地址的位开始输出。
	// {
	// 	printf("%d", (*p_data >> i) & 1);
	// }
	// printf("\n");
}

int to_absErrorBound_binary(double absErrBound)
{
  int n;
  for(n=0; n<100; n++)
  {
    if(absErrBound >= pow(2, -n))
    {
      return n;
    }
  }
}

uint32_t do_crc32(unsigned char *data_bits, int bytes)
{
    uint32_t crc = crc32(0L, Z_NULL, 0);

    for (int i = 0; i < bytes; i++)
    {
        crc = crc32(crc, data_bits + i, 1);
    }

    return crc;
}

uint64_t get_random_int(uint64_t from, uint64_t to)
{
  uint64_t n;
  n = rand() % (to - from + 1) + from;
  return n;
}
