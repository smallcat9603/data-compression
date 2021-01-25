#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//N is the max number of data bits
#define N 100

int hmLength(int k); //hamming check bits
void hamming_code(char* data, char* c, int k, int r); //value of each hamming check bit
void hamming_verify(char* data, char* c, int k, int r, char* v); //hamming verfify
int error_info(char* v, int r, int* error_bit_pos); //calculate bit error position if possible one bit error
void hamming_print(char* data, char* c, int k, int r); //print hamming secded
void hamming_rectify(char* data, char* c, int k, int r, int error_bit_pos); //rectify error bit if one error (except parity bit)
void cast_bits_to_char(unsigned char* bits, char* data, int bytes); //1 bit --> char '0' or '1' (8 bits)
void hamming_encode(unsigned char* bits, char** c, int bytes, int* r); //data bits --> char[]
int hamming_decode(unsigned char* bits, char* c, int bytes, int r);
void hamming_verify_bit(unsigned char* bits, char* c, int bytes, int r, char* v);
void hamming_rectify_bit(unsigned char* bits, char* c, int bytes, int r, int error_bit_pos); //rectify error bit if one error (except parity bit)

void main()
{
  int k = 0, r = 0, dnum = 0, cnum = 0; //1100 -> 0111100, k = 4, r = 3
  char data[N];
  char c[N];

  //test
  // data[0] = 'z';
  // data[1] = '7';
  // char* d = NULL;
  // hamming_encode(data, &d, 2, &r);
  // for(int i=0; i<(r+1); i++)
  // {
  //   printf("%c", d[i]);
  // }
  // printf("\n");
  // data[0] ^= 1 << 7;
  // c[0] = '0';
  // c[1] = '1';
  // c[2] = '0';
  // c[3] = '1';
  // c[4] = '0';
  // c[5] = '0';
  // hamming_decode(data, c, 2, 5);
  // return;

  printf("Input data to encode: ");

  for(k = 0; k < N; k++)
  {
    data[k] = getchar();
    if(data[k] != '0' && data[k] != '1') break;
  }

  r = hmLength(k);
  hamming_code(data, c, k, r);
  printf("Hamming code is ");
  for(int j = 1; j < r+k+1; j++)
  {
    if(j == (int)pow(2, cnum))
    {
      printf("%c", c[cnum]);
      cnum++;
    }
    else
    {
      printf("%c", data[dnum]);
      dnum++;
    }
  }
  printf(" %c", c[r]);
  printf(" (%d bits, k = %d, r = %d)\n", k+r, k, r);

  // test
  // data[2] = data[2] == '0'?'1':'0';
  // data[3] = data[3] == '0'?'1':'0';
  // c[r-1] = c[r-1] == '0'?'1':'0';
  // c[r] = c[r] == '0'?'1':'0';

  printf("Hamming recv is ");
  hamming_print(data, c, k, r);

  //hamming verify
  char v[r+1];
  hamming_verify(data, c, k, r, v);
  printf("Hamming veri is ");
  for(int i = 0; i < r+1; i++)
  {
    if(i == r) printf(" ");
    printf("%c", v[i]);
  }
  int error_bit_pos = 0; // 1 for most left bit
  int error_type = error_info(v, r, &error_bit_pos);
  switch(error_type)
  {
    case 0:
      printf(" (no error)\n");
      break;
    case 1:
      printf(" (two-bit error)\n");
      break;
    case 2:
      printf(" (parity error)\n");
      c[r] = c[r] == '0'?'1':'0';
      break;
    case 3:
      printf(" (one bit error: pos = %d)\n", error_bit_pos);
      hamming_rectify(data, c, k, r, error_bit_pos);
      break;
    default:
      printf(" ERROR");
  }
  printf("Hamming rect is ");
  hamming_print(data, c, k, r);
}

//length(c) = r+1 (secded)
void hamming_code(char* data, char* c, int k, int r)
{
  for(int i=0; i<r; i++)
  {
    int sum = 0, dnum = 0, cnum = 0;
    for(int j=1; j<r+k+1; j++)
    {
      if(j == (int)pow(2, cnum))
      {
        cnum++;
      }
      else
      {
        int x = pow(2, i);
        int y = j%(x*2);
        x = y/x;
        sum += (data[dnum]-'0')*x;
        dnum++;
      }
    }
    c[i] = sum%2 == 0?'0':'1';
  }

  //one additional parity bit for two bit error detection (secded)
  int sum = 0;
  for(int i=0; i<k; i++)
  {
    sum += data[i] - '0';
  }
  for(int i=0; i<r; i++)
  {
    sum += c[i] - '0';
  }
  c[r] = sum%2 + '0';
}

//calculate r
int hmLength(int k)
{
  int r = 0, flag = 1;
  while(flag)
  {
    int temp = pow(2, r);
    temp = temp - 1;
    flag = (temp-r-k<0);
    r++;
  }
  return r-1;
}

//hamming verify --> v
void hamming_verify(char* data, char* c, int k, int r, char* v)
{
  for(int i=0; i<r; i++)
  {
    int sum = 0, dnum = 0, cnum = 0;
    for(int j=1; j<r+k+1; j++)
    {
      if(j == (int)pow(2, cnum))
      {
        cnum++;
      }
      else
      {
        int x = pow(2, i);
        int y = j%(x*2);
        x = y/x;
        sum += (data[dnum]-'0')*x;
        dnum++;
      }
    }
    v[i] = sum%2 == (c[i]-'0')?'0':'1';
  }

  int sum = 0;
  for(int i=0; i<k; i++)
  {
    sum += data[i] - '0';
  }
  for(int i=0; i<r; i++)
  {
    sum += c[i] - '0';
  }
  v[r] = sum%2 == (c[r]-'0')?'0':'1';
}

//calculate bit error position if possible one bit error
int error_info(char* v, int r, int* error_bit_pos)
{
  int error_type = 0;

  for(int i=0; i<r; i++)
  {
    *error_bit_pos += (v[i]-'0')*pow(2, i);
  }

  if(*error_bit_pos > 0 && v[r] == '0') //two bit error
  {
    error_type = 1;
  }
  else if(*error_bit_pos == 0 && v[r] == '1') //parity bit error
  {
    error_type = 2;
  }
  else if(*error_bit_pos > 0 && v[r] == '1') //one bit error
  {
    error_type = 3;
  }

  return error_type;
}

//print hamming secded
void hamming_print(char* data, char* c, int k, int r)
{
  int dnum = 0, cnum = 0; 

  for(int j = 1; j < r+k+1; j++)
  {
    if(j == (int)pow(2, cnum))
    {
      printf("%c", c[cnum]);
      cnum++;
    }
    else
    {
      printf("%c", data[dnum]);
      dnum++;
    }
  }
  printf(" %c\n", c[r]);  
}

//rectify error bit if one error (except parity bit)
void hamming_rectify(char* data, char* c, int k, int r, int error_bit_pos)
{
  int dnum = 0, cnum = 0;

  for(int j=1; j<r+k+1; j++)
  {
    if(j == (int)pow(2, cnum))
    {   
      if(j == error_bit_pos)
      {
        c[cnum] = c[cnum] == '0'?'1':'0';
        break;
      }
      else
      {
        cnum++;
      }  
    }
    else
    {  
      if(j == error_bit_pos)
      {
        data[dnum] = data[dnum] == '0'?'1':'0';
        break;
      }
      else
      {
        dnum++;
      } 
    }
  }
}

//1 bit --> char '0' or '1' (8 bits)
void cast_bits_to_char(unsigned char* bits, char* data, int bytes)
{
  for(int i = 0; i < bytes; i++)
  {
    for(int j = 0; j < 8; j++)
    {
      int pos = i * 8 + j;
      data[pos] = ((bits[i] >> (7-j)) & 1) + '0';
    }
  }

  // test
  for(int i = 0; i < bytes; i++)
  {
    for(int j = 0; j < 8; j++)
    {
      printf("%c", ((bits[i] >> (7-j)) & 1) + '0');
    }
  }
  printf("\n");
  for(int i = 0; i < bytes * 8; i++)
  {
    printf("%c", data[i]);
  }
  printf("\n");
}

//data bits --> char[]
void hamming_encode(unsigned char* bits, char** c, int bytes, int* r)
{
  char data[bytes*8];
  cast_bits_to_char(bits, data, bytes);

  *r = hmLength(bytes*8);
  *c = (char*) malloc(sizeof(char)*(*r+1));
  hamming_code(data, *c, bytes*8, *r);
}

int hamming_decode(unsigned char* bits, char* c, int bytes, int r)
{
  char v[r+1];
  hamming_verify_bit(bits, c, bytes, r, v);

  int error_bit_pos = 0; // 1 for most left bit
  int error_type = error_info(v, r, &error_bit_pos);  
  switch(error_type)
  {
    case 0:
      printf("no error\n");
      break;
    case 1:
      printf("two-bit error\n");
      break;
    case 2:
      printf("parity error\n");
      c[r] = c[r] == '0'?'1':'0';
      break;
    case 3:
      printf("one bit error: pos = %d\n", error_bit_pos);
      hamming_rectify_bit(bits, c, bytes, r, error_bit_pos);
      break;
    default:
      printf("ERROR");
  }

  return error_type;
}

//hamming verify --> v
void hamming_verify_bit(unsigned char* bits, char* c, int bytes, int r, char* v)
{
  for(int i=0; i<r; i++)
  {
    int sum = 0, bnum = 0, cnum = 0;
    for(int j=1; j<r+bytes*8+1; j++)
    {
      if(j == (int)pow(2, cnum))
      {
        cnum++;
      }
      else
      {
        int x = pow(2, i);
        int y = j%(x*2);
        x = y/x;
        int bits_num = bnum/8;
        int bit_num = bnum%8;       
        sum += ((bits[bits_num] >> (7 - bit_num)) & 1) * x;
        bnum++;
      }
    }
    v[i] = sum%2 == (c[i]-'0')?'0':'1';
  }

  int sum = 0;
  for(int i=0; i<bytes; i++)
  {
    for(int j=0; j<8; j++)
    {
      sum += (bits[i] >> j) & 1;
    }
  }
  for(int i=0; i<r; i++)
  {
    sum += c[i] - '0';
  }
  v[r] = sum%2 == (c[r]-'0')?'0':'1';
}

//rectify error bit if one error (except parity bit)
void hamming_rectify_bit(unsigned char* bits, char* c, int bytes, int r, int error_bit_pos)
{
  int bnum = 0, cnum = 0;

  for(int j=1; j<r+bytes*8+1; j++)
  {
    if(j == (int)pow(2, cnum))
    {   
      if(j == error_bit_pos)
      {
        c[cnum] = c[cnum] == '0'?'1':'0';
        break;
      }
      else
      {
        cnum++;
      }  
    }
    else
    {  
      if(j == error_bit_pos)
      {
        int bits_num = bnum/8;
        int bit_num = bnum%8;
        bits[bits_num] ^= 1 << (7 - bit_num);
        break;
      }
      else
      {
        bnum++;
      } 
    }
  }
}