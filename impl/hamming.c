#include <stdio.h>
#include <math.h>

//N is the max number of data bits
#define N 100

int hmLength(int k); //hamming check bits
void hamming_code(char* data, char* c, int k, int r); //value of each hamming check bit

void main()
{
  int k = 0, r = 0, dnum = 0, cnum = 0; //1100 -> 0111100, k = 4, r = 3
  char data[N];
  char c[N];

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
  printf(" (%d bits, k = %d, r = %d)\n", k+r, k, r);
}

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
        sum += data[dnum]*x;
        dnum++;
      }
    }
    c[i] = sum%2 == 0?'0':'1';
  }
}

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
