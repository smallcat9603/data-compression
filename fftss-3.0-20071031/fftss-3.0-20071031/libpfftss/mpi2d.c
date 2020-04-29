/*
 * Copyright (C) 2002-2007 Akira Nukada. All rights reserved.
 * Copyright (C) 2002-2007 The SSI Project. All rights reserved. 
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the project nor the names of its contributors 
 *       may be used to endorse or promote products derived from this software 
 *       without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE SSI PROJECT ``AS IS'' AND ANY EXPRESS 
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE SSI PROJECT BE LIABLE FOR ANY DIRECT,
 * INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
 *  STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#elif defined(_WIN32)
#include "win32config.h"
#endif
#include <stdio.h>
#include <malloc.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <mpi.h>
#include <fftss.h>

#ifndef _OPENMP
#define omp_get_max_threads()   1
#define omp_get_thread_num()    0
#endif

//todo
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#define absErrBound         0.000001 //default 0.0001=2^{-12} (-13?), 0.000001=2^{-20}, 0.00001=2^{-16}, 0.001=2^{-10}, 0.01=2^{-7}
#define absErrBound_binary  20 //bitwise, SZ, equal to above
#define CT                  5 //compress type for pingpong & himeno & k-means, 0 no compress, 1 mycompress, 2 no-lossy-performance, 3 no-lossy-area, 4 sz, 5 bitwise
#define byte_or_bit         2 //1 byte, 2 bit

int myCompress_double(double[], double**, char**, int**, int);
double* myDecompress_double(double[], char[], int[], int);
// static void alltoalls0c(pfftss_plan);
// static void alltoalls1c(pfftss_plan);
double* myDecompress_bitwise_double(unsigned char*, int, int);
double decompress_bitwise_double(char*, int, double, double, double);
void myCompress_bitwise_double(double[], int, unsigned char**, int*, int*);
void compress_bitwise_double(double, unsigned char**, int*, int*);
double toSmallDataset_double(double[], double**, int);
void doubletostr(double*, char*);
double strtodbl(char*);
void add_bit_to_bytes(unsigned char**, int*, int*, int);
void bit_set(unsigned char*, unsigned char, int);
float calcCompressionRatio_sz_double(double[], int);
float calcCompressionRatio_nolossy_performance_double(double[], int);
float calcCompressionRatio_nolossy_area_double(double[], int);
void getDoubleBin(double,char[]);

typedef struct _pfftss_plan_s {
  void (*fp)(struct _pfftss_plan_s *, double *);
  long nx, ny, oy, ly, ox, lx;
  int id, npe;
  long py, plx;
  double **b;
  double *inout, *sb, *rb;
  long *oxp, *oyp;
  int *soff, *roff;
  int *rsoff, *rroff;
  int *ssize, *rsize;
  MPI_Request *sreq, *rreq;
  MPI_Status *st;
  fftss_plan p1, p2;
  MPI_Comm comm;
#ifdef HAVE_MPI2
  MPI_Win rwin, swin;
#endif
  long flags;
} pfftss_plan_s;

typedef pfftss_plan_s *pfftss_plan;


#if 1
extern void fftss_copyin(double *, double *, long, long);
extern void fftss_copyout(double *, double *, long, long);
#endif

static void trans_x(pfftss_plan p, double *inout)
{
  long i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->ly; i++)
    fftss_execute_dft(p->p1, inout + i * p->py * 2, 
		      p->b[2 * omp_get_thread_num()]);
}

static void trans_y(pfftss_plan p)
{
  long i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->lx; i++) {
    double *b0, *b1;

    b0 = p->b[2 * omp_get_thread_num()];
    b1 = p->b[2 * omp_get_thread_num() + 1];
    fftss_copyin(b0, p->rb + i * 2, p->plx, p->ny);
    fftss_execute_dft(p->p2, b0, b1);
    fftss_copyout(p->rb + i * 2, b0, p->plx, p->ny);
  }
}

static void alltoall0(pfftss_plan p)
{
  MPI_Alltoall(p->sb, p->plx * p->ly * 2, MPI_DOUBLE,
	       p->rb, p->plx * p->ly * 2, MPI_DOUBLE,
	       p->comm);
}

static void alltoalls0(pfftss_plan p)
{
  long i;
  int bsize;

  bsize = p->plx * p->ly * 2;

  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(p->rb + bsize * d, bsize, MPI_DOUBLE, d, 
    MPI_ANY_TAG, p->comm, &(p->rreq[i - 1]));   
  }
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Isend(p->sb + bsize * d, bsize, MPI_DOUBLE, d, 
        0, p->comm, &(p->sreq[i - 1]));  
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->rb[bsize * p->id + i] = p->sb[bsize * p->id + i];

  MPI_Waitall(p->npe - 1, p->rreq, p->st);
  MPI_Waitall(p->npe - 1, p->sreq, p->st);    
}

static void alltoalls0c(pfftss_plan p)
{
  long i;
  int bsize;

  double time_compress = 0;
  double start_time;
  double end_time;
  double start_time0 = fftss_get_wtime();
  double end_time0;

  bsize = p->plx * p->ly * 2;

  float compress_ratio = 0;
  int array_double_len_send[p->npe - 1];
  int array_double_len_recv[p->npe - 1];

  double* array_double_send[p->npe - 1];
  char* array_char_send[p->npe - 1];
  int* array_char_displacement_send[p->npe - 1];

  MPI_Request* mr0 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms0 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(&array_double_len_recv[i - 1], 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2]);
    double send_data[bsize];
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      send_data[j] = p->sb[first_index++];
    }
    array_double_send[i - 1] = NULL;
    array_char_send[i - 1] = NULL;
    array_char_displacement_send[i - 1] = NULL;
    int array_double_len = myCompress_double(send_data, &array_double_send[i - 1], &array_char_send[i - 1], &array_char_displacement_send[i - 1], bsize);
    array_double_len_send[i - 1] = array_double_len;
    MPI_Isend(&array_double_len, 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2 + 1]); 
  }
  MPI_Waitall((p->npe - 1) * 2, mr0, ms0);
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;

  double* array_double_recv[p->npe - 1]; 
  char* array_char_recv[p->npe - 1]; 
  int* array_char_displacement_recv[p->npe - 1]; 

  MPI_Request* mr1 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 3); 
  MPI_Status* ms1 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 3); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;  
    int num_p = array_double_len_recv[i - 1], num_c = bsize - array_double_len_recv[i - 1];
    array_double_recv[i - 1] = (double*) malloc(sizeof(double)*num_p);
    array_char_recv[i - 1] = (char*) malloc(sizeof(char)*num_c);
    array_char_displacement_recv[i - 1] = (int*) malloc(sizeof(int)*num_c);
    MPI_Irecv(array_double_recv[i - 1], num_p, MPI_DOUBLE, d, 1, p->comm, &mr1[(i - 1) * 3]);
    MPI_Irecv(array_char_recv[i - 1], num_c, MPI_CHAR, d, 2, p->comm, &mr1[(i - 1) * 3 + 1]);
    MPI_Irecv(array_char_displacement_recv[i - 1], num_c, MPI_INT, d, 3, p->comm, &mr1[(i - 1) * 3 + 2]);    
  }

  MPI_Request* mr2 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 3); 
  MPI_Status* ms2 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 3); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    //mycompress
    int num_p = array_double_len_send[i - 1], num_c = bsize - array_double_len_send[i - 1];
    compress_ratio += (float)(num_c*sizeof(char)+num_p*sizeof(double))/((num_c+num_p)*sizeof(double)); 
    MPI_Isend(array_double_send[i - 1], num_p, MPI_DOUBLE, d, 1, p->comm, &mr2[(i - 1) * 3]); 
    MPI_Isend(array_char_send[i - 1], num_c, MPI_CHAR, d, 2, p->comm, &mr2[(i - 1) * 3 + 1]); 
    MPI_Isend(array_char_displacement_send[i - 1], num_c, MPI_INT, d, 3, p->comm, &mr2[(i - 1) * 3 + 2]); 
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->rb[bsize * p->id + i] = p->sb[bsize * p->id + i];

  MPI_Waitall((p->npe - 1) * 3, mr1, ms1);
  MPI_Waitall((p->npe - 1) * 3, mr2, ms2);  

  //mydecompress
  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    double* decompressed_data = myDecompress_double(array_double_recv[i - 1], array_char_recv[i - 1], array_char_displacement_recv[i - 1], bsize);
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      p->rb[first_index++] = decompressed_data[j];
    }
  }
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;
  end_time0 = fftss_get_wtime();

  if(p->id == 0) 
  {
    printf("compress ratio = %f \n", 1/(compress_ratio/(p->npe - 1)));
    printf("compress time, total time  = %f, %f \n", time_compress, end_time0 - start_time0);
  }
}

static void alltoalls0cb(pfftss_plan p)
{
  long i;
  int bsize;

  double time_compress = 0;
  double start_time;
  double end_time;
  double start_time0 = fftss_get_wtime();;
  double end_time0;

  bsize = p->plx * p->ly * 2;

  float compress_ratio = 0;

  int data_bytes_send[p->npe - 1];
  int data_bytes_recv[p->npe - 1];

  unsigned char* data_bits_send[p->npe - 1];
  double data_min_send[p->npe - 1];

  MPI_Request* mr0 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms0 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(&data_bytes_recv[i - 1], 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2]);
    double send_data[bsize];
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      send_data[j] = p->sb[first_index++];
    }

    float sz_comp_ratio = calcCompressionRatio_sz_double(send_data, bsize);
    float nolossy_performance = calcCompressionRatio_nolossy_performance_double(send_data, bsize);
    float nolossy_area = calcCompressionRatio_nolossy_area_double(send_data, bsize);
    printf("compression ratio: sz %f, nolossy_performance %f, nolossy_area %f \n", 1/sz_comp_ratio, 1/nolossy_performance, 1/nolossy_area);

    double* send_data_small = NULL;
    data_min_send[i - 1] = toSmallDataset_double(send_data, &send_data_small, bsize);          

    data_bits_send[i - 1] = NULL;
    int pos = 8; //position of filled bit in last byte --> 87654321

    myCompress_bitwise_double(send_data_small, bsize, &data_bits_send[i - 1], &data_bytes_send[i - 1], &pos);
    compress_ratio += data_bytes_send[i - 1]*8.0/(bsize*sizeof(double)*8); 
    MPI_Isend(&data_bytes_send[i - 1], 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2 + 1]); 
  }
  MPI_Waitall((p->npe - 1) * 2, mr0, ms0);
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;

  unsigned char* data_bits_recv[p->npe - 1];
  double data_min_recv[p->npe - 1];

  MPI_Request* mr1 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms1 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    data_bits_recv[i - 1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[i - 1]);
    MPI_Irecv(&data_min_recv[i - 1], 1, MPI_DOUBLE, d, 1, p->comm, &mr1[(i - 1) * 2]);
    MPI_Irecv(data_bits_recv[i - 1], data_bytes_recv[i - 1], MPI_UNSIGNED_CHAR, d, 2, p->comm, &mr1[(i - 1) * 2 + 1]);
  }

  MPI_Request* mr2 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms2 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    //bitwise compress
    MPI_Isend(&data_min_send[i - 1], 1, MPI_DOUBLE, d, 1, p->comm, &mr2[(i - 1) * 2]); 
    MPI_Isend(data_bits_send[i - 1], data_bytes_send[i - 1], MPI_UNSIGNED_CHAR, d, 2, p->comm, &mr2[(i - 1) * 2 + 1]); 
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->rb[bsize * p->id + i] = p->sb[bsize * p->id + i];

  MPI_Waitall((p->npe - 1) * 2, mr1, ms1);
  MPI_Waitall((p->npe - 1) * 2, mr2, ms2);  

  //bitwise decompress
  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    double* decompressed_data = myDecompress_bitwise_double(data_bits_recv[i - 1], data_bytes_recv[i - 1], bsize);
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      p->rb[first_index++] = decompressed_data[j] + data_min_recv[i - 1];
    }
  }
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;
  end_time0 = fftss_get_wtime();

  if(p->id == 0) 
  {
    printf("compress ratio (0cb) = %f \n", 1/(compress_ratio/(p->npe - 1)));
    printf("compress time, total time  = %f, %f \n", time_compress, end_time0 - start_time0);
  }
}

#ifdef HAVE_MPI2

static void alltoallp0(pfftss_plan p)
{
  long i;
  int bsize;

  MPI_Win_fence(MPI_MODE_NOPRECEDE, p->rwin);

  bsize = p->plx * p->ly * 2;

  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Put(p->sb + bsize * d, bsize, MPI_DOUBLE, d, 
	    bsize * p->id, bsize, MPI_DOUBLE, p->rwin);
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->rb[bsize * p->id + i] = p->sb[bsize * p->id + i];

  MPI_Win_fence(MPI_MODE_NOSUCCEED, p->rwin);
}

#endif

static void alltoallv0(pfftss_plan p)
{
  MPI_Alltoallv(p->sb, p->ssize, p->soff, MPI_DOUBLE,
		p->rb, p->rsize, p->roff, MPI_DOUBLE,
		p->comm);
}

static void alltoallvs0(pfftss_plan p)
{
  long i;

  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(p->rb + p->roff[d], p->rsize[d], MPI_DOUBLE, d, 
	      MPI_ANY_TAG, p->comm, &(p->rreq[i - 1]));
  }
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Isend(p->sb + p->soff[d], p->ssize[d], MPI_DOUBLE, d, 
	      0, p->comm, &(p->sreq[i - 1]));
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->ssize[p->id]; i++) 
    p->rb[p->roff[p->id] + i] = p->sb[p->soff[p->id] + i];

  MPI_Waitall(p->npe - 1, p->rreq, p->st);
  MPI_Waitall(p->npe - 1, p->sreq, p->st);
}

#ifdef HAVE_MPI2

static void alltoallvp0(pfftss_plan p)
{
  long i;

  MPI_Win_fence(MPI_MODE_NOPRECEDE, p->rwin);
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Put(p->sb + p->soff[d], p->ssize[d], MPI_DOUBLE, d,
	    p->rroff[d], p->ssize[d], MPI_DOUBLE, p->rwin);
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->ssize[p->id]; i++) 
    p->rb[p->roff[p->id] + i] = p->sb[p->soff[p->id] + i];

  MPI_Win_fence(MPI_MODE_NOSUCCEED, p->rwin);
}

#endif

static void alltoall1(pfftss_plan p)
{
  MPI_Alltoall(p->rb, p->plx * p->ly * 2, MPI_DOUBLE,
	       p->sb, p->plx * p->ly * 2, MPI_DOUBLE,
	       p->comm);
}

static void alltoalls1(pfftss_plan p)
{
  long i;
  int bsize;

  bsize = p->plx * p->ly * 2;

  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(p->sb + bsize * d, bsize, MPI_DOUBLE, d, 
        MPI_ANY_TAG, p->comm, &(p->rreq[i - 1]));      
  }
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Isend(p->rb + bsize * d, bsize, MPI_DOUBLE, d, 
        0, p->comm, &(p->sreq[i - 1]));      
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->sb[bsize * p->id + i] = p->rb[bsize * p->id + i];

  MPI_Waitall(p->npe - 1, p->rreq, p->st);
  MPI_Waitall(p->npe - 1, p->sreq, p->st);    
}

static void alltoalls1c(pfftss_plan p)
{
  long i;
  int bsize;

  double time_compress = 0;
  double start_time;
  double end_time;
  double start_time0 = fftss_get_wtime();
  double end_time0;

  bsize = p->plx * p->ly * 2;

  float compress_ratio = 0;
  int array_double_len_send[p->npe - 1];
  int array_double_len_recv[p->npe - 1];

  double* array_double_send[p->npe - 1];
  char* array_char_send[p->npe - 1];
  int* array_char_displacement_send[p->npe - 1];

  MPI_Request* mr0 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms0 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(&array_double_len_recv[i - 1], 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2]);
    double send_data[bsize];
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      send_data[j] = p->rb[first_index++];
    }
    array_double_send[i - 1] = NULL;
    array_char_send[i - 1] = NULL;
    array_char_displacement_send[i - 1] = NULL;
    int array_double_len = myCompress_double(send_data, &array_double_send[i - 1], &array_char_send[i - 1], &array_char_displacement_send[i - 1], bsize);
    array_double_len_send[i - 1] = array_double_len;
    MPI_Isend(&array_double_len, 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2 + 1]); 
  }
  MPI_Waitall((p->npe - 1) * 2, mr0, ms0);
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;

  double* array_double_recv[p->npe - 1]; 
  char* array_char_recv[p->npe - 1]; 
  int* array_char_displacement_recv[p->npe - 1]; 

  MPI_Request* mr1 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 3); 
  MPI_Status* ms1 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 3); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;  
    int num_p = array_double_len_recv[i - 1], num_c = bsize - array_double_len_recv[i - 1];
    array_double_recv[i - 1] = (double*) malloc(sizeof(double)*num_p);
    array_char_recv[i - 1] = (char*) malloc(sizeof(char)*num_c);
    array_char_displacement_recv[i - 1] = (int*) malloc(sizeof(int)*num_c);
    MPI_Irecv(array_double_recv[i - 1], num_p, MPI_DOUBLE, d, 1, p->comm, &mr1[(i - 1) * 3]);
    MPI_Irecv(array_char_recv[i - 1], num_c, MPI_CHAR, d, 2, p->comm, &mr1[(i - 1) * 3 + 1]);
    MPI_Irecv(array_char_displacement_recv[i - 1], num_c, MPI_INT, d, 3, p->comm, &mr1[(i - 1) * 3 + 2]);    
  }

  MPI_Request* mr2 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 3); 
  MPI_Status* ms2 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 3); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    //mycompress
    int num_p = array_double_len_send[i - 1], num_c = bsize - array_double_len_send[i - 1];
    compress_ratio += (float)(num_c*sizeof(char)+num_p*sizeof(double))/((num_c+num_p)*sizeof(double)); 
    MPI_Isend(array_double_send[i - 1], num_p, MPI_DOUBLE, d, 1, p->comm, &mr2[(i - 1) * 3]); 
    MPI_Isend(array_char_send[i - 1], num_c, MPI_CHAR, d, 2, p->comm, &mr2[(i - 1) * 3 + 1]); 
    MPI_Isend(array_char_displacement_send[i - 1], num_c, MPI_INT, d, 3, p->comm, &mr2[(i - 1) * 3 + 2]); 
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->sb[bsize * p->id + i] = p->rb[bsize * p->id + i];

  MPI_Waitall((p->npe - 1) * 3, mr1, ms1);
  MPI_Waitall((p->npe - 1) * 3, mr2, ms2); 

  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    double* decompressed_data = myDecompress_double(array_double_recv[i - 1], array_char_recv[i - 1], array_char_displacement_recv[i - 1], bsize);
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      p->sb[first_index++] = decompressed_data[j];
    }
  }
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;
  end_time0 = fftss_get_wtime();

  if(p->id == 0) 
  {
    printf("compress ratio = %f \n", 1/(compress_ratio/(p->npe - 1)));
    printf("compress time, total time  = %f, %f \n", time_compress, end_time0 - start_time0);
  }  
}

static void alltoalls1cb(pfftss_plan p)
{
  long i;
  int bsize;

  double time_compress = 0;
  double start_time;
  double end_time;
  double start_time0 = fftss_get_wtime();;
  double end_time0;

  bsize = p->plx * p->ly * 2;

  float compress_ratio = 0;

  int data_bytes_send[p->npe - 1];
  int data_bytes_recv[p->npe - 1];

  unsigned char* data_bits_send[p->npe - 1];
  double data_min_send[p->npe - 1];

  MPI_Request* mr0 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms0 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(&data_bytes_recv[i - 1], 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2]);
    double send_data[bsize];
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      send_data[j] = p->rb[first_index++];
    }

    float sz_comp_ratio = calcCompressionRatio_sz_double(send_data, bsize);
    float nolossy_performance = calcCompressionRatio_nolossy_performance_double(send_data, bsize);
    float nolossy_area = calcCompressionRatio_nolossy_area_double(send_data, bsize);
    printf("compression ratio: sz %f, nolossy_performance %f, nolossy_area %f \n", 1/sz_comp_ratio, 1/nolossy_performance, 1/nolossy_area);

    double* send_data_small = NULL;
    data_min_send[i - 1] = toSmallDataset_double(send_data, &send_data_small, bsize);    

    data_bits_send[i - 1] = NULL;
    int pos = 8; //position of filled bit in last byte --> 87654321

    myCompress_bitwise_double(send_data_small, bsize, &data_bits_send[i - 1], &data_bytes_send[i - 1], &pos);
    compress_ratio += data_bytes_send[i - 1]*8.0/(bsize*sizeof(double)*8); 
    MPI_Isend(&data_bytes_send[i - 1], 1, MPI_INT, d, 0, p->comm, &mr0[(i - 1) * 2 + 1]); 
  }
  MPI_Waitall((p->npe - 1) * 2, mr0, ms0);
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;

  unsigned char* data_bits_recv[p->npe - 1];
  double data_min_recv[p->npe - 1];

  MPI_Request* mr1 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms1 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    data_bits_recv[i - 1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[i - 1]);
    MPI_Irecv(&data_min_recv[i - 1], 1, MPI_DOUBLE, d, 1, p->comm, &mr1[(i - 1) * 2]);
    MPI_Irecv(data_bits_recv[i - 1], data_bytes_recv[i - 1], MPI_UNSIGNED_CHAR, d, 2, p->comm, &mr1[(i - 1) * 2 + 1]);
  }

  MPI_Request* mr2 = (MPI_Request *)malloc(sizeof(MPI_Request) * (p->npe - 1) * 2); 
  MPI_Status* ms2 = (MPI_Status *)malloc(sizeof(MPI_Status) * (p->npe - 1) * 2); 
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    //bitwise compress
    MPI_Isend(&data_min_send[i - 1], 1, MPI_DOUBLE, d, 1, p->comm, &mr2[(i - 1) * 2]); 
    MPI_Isend(data_bits_send[i - 1], data_bytes_send[i - 1], MPI_UNSIGNED_CHAR, d, 2, p->comm, &mr2[(i - 1) * 2 + 1]); 
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->sb[bsize * p->id + i] = p->rb[bsize * p->id + i];

  MPI_Waitall((p->npe - 1) * 2, mr1, ms1);
  MPI_Waitall((p->npe - 1) * 2, mr2, ms2);  

  //bitwise decompress
  start_time = fftss_get_wtime();
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    double* decompressed_data = myDecompress_bitwise_double(data_bits_recv[i - 1], data_bytes_recv[i - 1], bsize);
    int first_index = bsize * d;
    for(int j = 0; j < bsize; j++)
    {
      p->sb[first_index++] = decompressed_data[j] + data_min_recv[i - 1];
    }
  }
  end_time = fftss_get_wtime();
  time_compress += end_time - start_time;
  end_time0 = fftss_get_wtime();

  if(p->id == 0) 
  {
    printf("compress ratio (1cb) = %f \n", 1/(compress_ratio/(p->npe - 1)));
    printf("compress time, total time  = %f, %f \n", time_compress, end_time0 - start_time0);
  }
}

#ifdef HAVE_MPI2

static void alltoallp1(pfftss_plan p)
{
  long i;
  int bsize;

  MPI_Win_fence(MPI_MODE_NOPRECEDE, p->swin);

  bsize = p->plx * p->ly * 2;

  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Put(p->rb + bsize * d, bsize, MPI_DOUBLE, d, 
	    bsize * p->id, bsize, MPI_DOUBLE, p->swin);
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < bsize; i++) 
    p->sb[bsize * p->id + i] = p->rb[bsize * p->id + i];

  MPI_Win_fence(MPI_MODE_NOSUCCEED, p->swin);
}

#endif

static void alltoallv1(pfftss_plan p)
{
  MPI_Alltoallv(p->rb, p->rsize, p->roff, MPI_DOUBLE,
		p->sb, p->ssize, p->soff, MPI_DOUBLE,
		p->comm);
}

static void alltoallvs1(pfftss_plan p)
{
  long i;

  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Irecv(p->sb + p->soff[d], p->ssize[d], MPI_DOUBLE, d, 
	      MPI_ANY_TAG, p->comm, &(p->rreq[i - 1]));
  }
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Isend(p->rb + p->roff[d], p->rsize[d], MPI_DOUBLE, d, 
	      0, p->comm, &(p->sreq[i - 1]));
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->rsize[p->id]; i++) 
    p->sb[p->soff[p->id] + i] = p->rb[p->roff[p->id] + i];

  MPI_Waitall(p->npe - 1, p->rreq, p->st);
  MPI_Waitall(p->npe - 1, p->sreq, p->st);
}

#ifdef HAVE_MPI2

static void alltoallvp1(pfftss_plan p)
{
  long i;

  MPI_Win_fence(MPI_MODE_NOPRECEDE, p->swin);
  for (i = 1; i < p->npe; i++) {
    int d;
    d = (p->id + i) % p->npe;
    MPI_Put(p->rb + p->roff[d], p->rsize[d], MPI_DOUBLE, d, 
	    p->rsoff[d], p->rsize[d], MPI_DOUBLE, p->swin);
  }

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->rsize[p->id]; i++) 
    p->sb[p->soff[p->id] + i] = p->rb[p->roff[p->id] + i];

  MPI_Win_fence(MPI_MODE_NOSUCCEED, p->swin);
}

#endif

static void pack0(pfftss_plan p, double *inout)
{
  long i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->npe; i++) {
    long j, k;
    long lxp;
    double *src, *dst;

    lxp = p->plx;
    for (j = 0; j < p->ly; j++) {
      src = inout + p->lx * i * 2 + p->py * j * 2;
      dst = p->sb + p->ly * lxp * i * 2 + lxp * j * 2;
      for (k = 0; k < p->lx * 2; k++)
	dst[k] = src[k];
    }
  }
}  

static void packv0(pfftss_plan p, double *inout)
{
  long i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->npe; i++) {
    long j, k;
    long lx, lxp;
    double *src, *dst;

    lx = p->oxp[i + 1] - p->oxp[i];
    lxp = p->plx;
    for (j = 0; j < p->ly; j++) {
      src = inout + p->oxp[i] * 2 + p->py * j * 2;
      dst = p->sb + p->soff[i] + lxp * j * 2;
      for (k = 0; k < lx * 2; k++)
	dst[k] = src[k];
    }
  }
}  

static void unpack1(pfftss_plan p, double *inout)
{
  long i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->npe; i++) {
    long j, k;
    long lxp;
    double *src, *dst;

    lxp = p->plx;
    for (j = 0; j < p->ly; j++) {
      dst = inout + p->lx * i * 2 + p->py * j * 2;
      src = p->sb + lxp * p->ly * i * 2 + lxp * j * 2;
      for (k = 0; k < p->lx * 2; k++)
	dst[k] = src[k];
    }
  }
}

static void unpackv1(pfftss_plan p, double *inout)
{
  long i;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for (i = 0; i < p->npe; i++) {
    long j, k;
    long lx, lxp;
    double *src, *dst;

    lx = p->oxp[i + 1] - p->oxp[i];
    lxp = p->plx;
    for (j = 0; j < p->ly; j++) {
      dst = inout + p->oxp[i] * 2 + p->py * j * 2;
      src = p->sb + p->soff[i] + lxp * j * 2;
      for (k = 0; k < lx * 2; k++)
	dst[k] = src[k];
    }
  }
}

void pfftss_execute_pdft_2d_a(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  pack0(p, inout);
  alltoall0(p);
  trans_y(p);
  alltoall1(p);
  unpack1(p, inout);

  return;
}

void pfftss_execute_pdft_2d_s(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  pack0(p, inout);
  alltoalls0(p);
  trans_y(p);
  alltoalls1(p);
  unpack1(p, inout);

  return;
}

void pfftss_execute_pdft_2d_sc(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  pack0(p, inout);
  alltoalls0c(p);
  trans_y(p);
  alltoalls1c(p);
  unpack1(p, inout);

  return;
}

void pfftss_execute_pdft_2d_scb(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  pack0(p, inout);
  alltoalls0cb(p);
  trans_y(p);
  alltoalls1cb(p);
  unpack1(p, inout);

  return;
}

#ifdef HAVE_MPI2

void pfftss_execute_pdft_2d_p(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  pack0(p, inout);
  alltoallp0(p);
  trans_y(p);
  alltoallp1(p);
  unpack1(p, inout);

  return;
}

#endif

void pfftss_execute_pdft_2d_va(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  packv0(p, inout);
  alltoallv0(p);
  trans_y(p);
  alltoallv1(p);
  unpackv1(p, inout);

  return;
}

void pfftss_execute_pdft_2d_vs(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  packv0(p, inout);
  alltoallvs0(p);
  trans_y(p);
  alltoallvs1(p);
  unpackv1(p, inout);

  return;
}

#ifdef HAVE_MPI2

void pfftss_execute_pdft_2d_vp(pfftss_plan p, double *inout)
{
  trans_x(p, inout);
  packv0(p, inout);
  alltoallvp0(p);
  trans_y(p);
  alltoallvp1(p);
  unpackv1(p, inout);

  return;
}

#endif

void pfftss_execute(pfftss_plan p)
{
  (*(p->fp))(p, p->inout);
  return;
}

void pfftss_execute_dft(pfftss_plan p, double *inout)
{
  (*(p->fp))(p, inout);
  return;
}

void pfftss_destroy_plan(pfftss_plan p)
{
  long i;

#ifdef HAVE_MPI_ALLOC_MEM
  MPI_Free_mem(p->sb);
  MPI_Free_mem(p->rb);
#else
  fftss_free(p->sb);
  fftss_free(p->rb);
#endif
  free(p->oxp);
  free(p->oyp);

  if (p->soff) {
    free(p->soff);
    free(p->roff);
    free(p->rsoff);
    free(p->rroff);
    free(p->ssize);
    free(p->rsize);
  }
  free(p->sreq);
  free(p->rreq);
  free(p->st);

  for (i = 0; i < omp_get_max_threads() * 2; i++)
    fftss_free(p->b[i]);
  free(p->b);
  if (p->p1 != p->p2) fftss_destroy_plan(p->p2);
  fftss_destroy_plan(p->p1);

#ifdef HAVE_MPI2
  MPI_Win_free(&p->rwin);
  MPI_Win_free(&p->swin);
#endif

  free(p);
}

static double test_run(pfftss_plan p, void (*fp)(pfftss_plan, double *))
{
  double t;

  MPI_Barrier(p->comm);
  t = fftss_get_wtime();
  (*fp)(p, p->inout);
  MPI_Barrier(p->comm);
  
  t = fftss_get_wtime() - t;
  if ((p->flags & FFTSS_VERBOSE) && p->id == 0) printf("%lf sec.\n", t);
  return t;
}

pfftss_plan 
pfftss_plan_dft_2d(long nx, long ny, long py, long oy, long ly,
		   double *inout, long sign, long flags, MPI_Comm comm)
{
  pfftss_plan p;
  long l;
  long i;
  double t0, t1;

  p = (pfftss_plan)malloc(sizeof(pfftss_plan_s));
  p->nx = nx; p->ny = ny; p->py = py; 
  p->oy = oy; p->ly = ly;
  p->comm = comm;
  p->flags = flags;
  MPI_Comm_size(comm, &(p->npe));
  MPI_Comm_rank(comm, &(p->id));

  p->ox = nx * p->id / p->npe;
  p->lx = nx * (p->id + 1) / p->npe - p->ox;
  p->plx = (nx + p->npe - 1) / p->npe;
  if ((p->plx & 1) == 0) p->plx++;

  p->inout = inout;

  l = (nx > ny) ? nx : ny;

  p->b = (double **)malloc(sizeof(double *) * 2 * omp_get_max_threads());

  for (i = 0; i < omp_get_max_threads(); i++) {
    p->b[i * 2] = (double *)fftss_malloc(sizeof(double) * 2 * l);
    p->b[i * 2 + 1] = (double *)fftss_malloc(sizeof(double) * 2 * l);
  }

  p->sreq = (MPI_Request *)malloc(sizeof(MPI_Request) * p->npe); 
  p->rreq = (MPI_Request *)malloc(sizeof(MPI_Request) * p->npe); 
  p->st = (MPI_Status *)malloc(sizeof(MPI_Status) * p->npe); 

#ifdef _OPENMP
#pragma omp parallel private(i)
#endif
  {
    double *b0, *b1;
    int id;

    id = omp_get_thread_num();
    b0 = p->b[id * 2];
    b1 = p->b[id * 2 + 1];
    for (i = 0; i < l * 2; i++) b0[i] = 0.0;
    for (i = 0; i < l * 2; i++) b1[i] = 0.0;
  } 
  
  p->p1 = fftss_plan_dft_1d(nx, p->b[0], p->b[1], sign, FFTSS_INOUT);
  if (nx != ny)
    p->p2 = fftss_plan_dft_1d(ny, p->b[0], p->b[1], sign, FFTSS_INOUT);
  else p->p2 = p->p1;

  p->oxp = (long *)malloc(sizeof(long) * (p->npe + 1));
  p->oyp = (long *)malloc(sizeof(long) * (p->npe + 1));

  MPI_Allgather(&(p->ox), 1, MPI_LONG, p->oxp, 1, MPI_LONG, comm);
  MPI_Allgather(&(p->oy), 1, MPI_LONG, p->oyp, 1, MPI_LONG, comm);

  p->oxp[p->npe] = nx;
  p->oyp[p->npe] = ny;

  l = p->ly;
  for (i = 0; i < p->npe; i++)
    if (p->oyp[i + 1] - p->oyp[i] != l) break;
  if (nx % p->npe == 0 && ny % p->npe == 0 && i == p->npe) {

    if ((p->flags & FFTSS_VERBOSE) && p->id == 0) printf("uniform mode\n");

    p->soff = NULL;
    p->roff = NULL;
    p->rsoff = NULL;
    p->rroff = NULL;
    p->ssize = NULL;
    p->rsize = NULL;

#ifdef HAVE_MPI_ALLOC_MEM
    MPI_Alloc_mem(p->plx * p->ny * sizeof(double) * 2,
		  MPI_INFO_NULL, &(p->sb));
    MPI_Alloc_mem(p->plx * p->ny * sizeof(double) * 2,
		  MPI_INFO_NULL, &(p->rb));
#else
    p->sb = 
      (double *)fftss_malloc(p->plx * p->ny * sizeof(double) * 2);
    p->rb =
      (double *)fftss_malloc(p->plx * p->ny * sizeof(double) * 2);
#endif

#ifdef HAVE_MPI2
    MPI_Win_create(p->sb, 
		   p->plx * p->ny * sizeof(double) * 2,
		   sizeof(double), MPI_INFO_NULL, comm, &(p->swin));
    MPI_Win_create(p->rb, 
		   p->plx * p->ny * sizeof(double) * 2,
		   sizeof(double), MPI_INFO_NULL, comm, &(p->rwin));
#endif

    if(CT == 0)
    {
      p->fp = pfftss_execute_pdft_2d_s;
      test_run(p, pfftss_execute_pdft_2d_s);      
    }
    else if(CT == 1)
    {
      p->fp = pfftss_execute_pdft_2d_sc;
      test_run(p, pfftss_execute_pdft_2d_sc);    
    }
    else if(CT == 5)
    {
      p->fp = pfftss_execute_pdft_2d_scb;
      test_run(p, pfftss_execute_pdft_2d_scb);          
    }
    else
    {
      p->fp = pfftss_execute_pdft_2d_a;
      t0 = test_run(p, pfftss_execute_pdft_2d_a);
      t1 = test_run(p, pfftss_execute_pdft_2d_a);
      if (t1 < t0) { t0 = t1; }
      t1 = test_run(p, pfftss_execute_pdft_2d_s);
      if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_s; t0 = t1; }
      t1 = test_run(p, pfftss_execute_pdft_2d_s);
      if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_s; t0 = t1; }    
    }
    

#ifdef HAVE_MPI2
    t1 = test_run(p, pfftss_execute_pdft_2d_p);
    if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_p; t0 = t1; }
    t1 = test_run(p, pfftss_execute_pdft_2d_p);
    if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_p; t0 = t1; }
#endif

  } else {

    p->ssize = (int *)malloc(sizeof(int) * p->npe);
    p->rsize = (int *)malloc(sizeof(int) * p->npe);
    p->soff = (int *)malloc(sizeof(int) * (p->npe + 1));
    p->roff = (int *)malloc(sizeof(int) * (p->npe + 1));
    
    for (i = 0; i < p->npe; i++) {
      p->ssize[i] = p->plx * ly * 2;
      p->rsize[i] = (p->oyp[i + 1] - p->oyp[i]) * p->plx * 2;
    }
    
    p->soff[0] = 0; p->roff[0] = 0;
    for (i = 0; i < p->npe; i++) {
      p->soff[i + 1] = p->soff[i] + p->ssize[i];
      p->roff[i + 1] = p->roff[i] + p->rsize[i];
    }

#ifdef HAVE_MPI_ALLOC_MEM    
    MPI_Alloc_mem(p->soff[p->npe] * sizeof(double), MPI_INFO_NULL, &(p->sb));
    MPI_Alloc_mem(p->roff[p->npe] * sizeof(double), MPI_INFO_NULL, &(p->rb));
#else
    p->sb = (double *)fftss_malloc(p->soff[p->npe] * sizeof(double));
    p->rb = (double *)fftss_malloc(p->roff[p->npe] * sizeof(double));
#endif

    p->rroff = (int *)malloc(sizeof(int) * p->npe);
    p->rsoff = (int *)malloc(sizeof(int) * p->npe);
    MPI_Alltoall(p->roff, 1, MPI_INT, p->rroff, 1, MPI_INT, comm);
    MPI_Alltoall(p->soff, 1, MPI_INT, p->rsoff, 1, MPI_INT, comm);

#ifdef HAVE_MPI2
    MPI_Win_create(p->sb, p->soff[p->npe] * sizeof(double), sizeof(double),
		   MPI_INFO_NULL, comm, &(p->swin));
    MPI_Win_create(p->rb, p->roff[p->npe] * sizeof(double), sizeof(double),
		   MPI_INFO_NULL, comm, &(p->rwin));
#endif
    
    p->fp = pfftss_execute_pdft_2d_va;
    t0 = test_run(p, pfftss_execute_pdft_2d_va);
    t1 = test_run(p, pfftss_execute_pdft_2d_va);
    if (t1 < t0) { t0 = t1; }
    t1 = test_run(p, pfftss_execute_pdft_2d_vs);
    if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_vs; t0 = t1; }
    t1 = test_run(p, pfftss_execute_pdft_2d_vs);
    if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_vs; t0 = t1; }

#ifdef HAVE_MPI2
    t1 = test_run(p, pfftss_execute_pdft_2d_vp);
    if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_vp; t0 = t1; }
    t1 = test_run(p, pfftss_execute_pdft_2d_vp);
    if (t1 < t0) { p->fp = pfftss_execute_pdft_2d_vp; t0 = t1; }
#endif
  }

  return p;
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
        printf("Error (re)allocating memory1");
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
      
      if(diff_min<=absErrBound) 
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
          printf("Error (re)allocating memory2");
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
          printf("Error (re)allocating memory3");
          exit(1);
        }             
      }
    }
  }
  return array_double_len;
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
              printf("Error (re)allocating memory4\n");
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
              printf("Error (re)allocating memory5\n");
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

          int mantissa_bits_within_error_bound = absErrBound_binary + expo_value;
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
              printf("Error (re)allocating memory6\n");
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
                printf("Error (re)allocating memory7\n");
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
                printf("Error (re)allocating memory8\n");
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
          printf("Error (re)allocating memory9\n");
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
      if(fabs(real_value) < absErrBound)
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
      
      if(fabs(real_value) < absErrBound)
      {
        add_bit_to_bytes(data_bits, bytes, pos, 1);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        add_bit_to_bytes(data_bits, bytes, pos, 0);
        d++;
      }
      else if(diff_min<=absErrBound) 
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

  int mantissa_bits_within_error_bound = absErrBound_binary + expo_value;

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
        mantissa_bits_within_error_bound = absErrBound_binary + expo_value;
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
        printf("Error (re)allocating memory0\n");
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

void getDoubleBin(double num,char bin[])
{
    int t = 1;
    int *f = (int*)(&num);
    for(int i=0;i<64;i++)
    {
        bin[i] = (*f)&(t<<63-i)?1:0;
    }
}
