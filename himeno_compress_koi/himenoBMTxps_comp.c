/********************************************************************
 
 modified by huyao for data compression

 This benchmark test program is measuring a cpu performance
 of floating point operation by a Poisson equation solver.

 If you have any question, please ask me via email.
 written by Ryutaro HIMENO, November 26, 2001.
 Version 3.0
 ----------------------------------------------
 Ryutaro Himeno, Dr. of Eng.
 Head of Computer Information Division,
 RIKEN (The Institute of Pysical and Chemical Research)
 Email : himeno@postman.riken.go.jp
 ---------------------------------------------------------------
 You can adjust the size of this benchmark code to fit your target
 computer. In that case, please chose following sets of
 (mimax,mjmax,mkmax):
 small : 33,33,65
 small : 65,65,129
 midium: 129,129,257
 large : 257,257,513
 ext.large: 513,513,1025
 This program is to measure a computer performance in MFLOPS
 by using a kernel which appears in a linear solver of pressure
 Poisson eq. which appears in an incompressible Navier-Stokes solver.
 A point-Jacobi method is employed in this solver as this method can 
 be easyly vectrized and be parallelized.
 ------------------
 Finite-difference method, curvilinear coodinate system
 Vectorizable and parallelizable on each grid point
 No. of grid points : imax x jmax x kmax including boundaries
 ------------------
 A,B,C:coefficient matrix, wrk1: source term of Poisson equation
 wrk2 : working area, OMEGA : relaxation parameter
 BND:control variable for boundaries and objects ( = 0 or 1)
 P: pressure
********************************************************************/

#include <stdio.h>
#include "mpi.h"
#include "param.h"
#include "dataCompression.h"
#include <stdlib.h>
#include <math.h>
//#include <stdint.h>
//#include <stdbool.h>
#include <assert.h>


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
static float  bnd[MIMAX][MJMAX][MKMAX];
static float  wrk1[MIMAX][MJMAX][MKMAX],
              wrk2[MIMAX][MJMAX][MKMAX];
static float omega;
static int npe,id;

static int ndx,ndy,ndz;
static int imax,jmax,kmax;
static int ndims=3,iop[3];
static int npx[2],npy[2],npz[2];
MPI_Comm     mpi_comm_cart;
MPI_Datatype ijvec,ikvec,jkvec;

//todo
static float cr = 0; //compression rate
static int cr_num = 0;

int
main(int argc,char *argv[])
{
  int    i,j,k,nn;
  int    mx,my,mz,it;
  float  gosa;
  double cpu,cpu0,cpu1,flop,target;
  target= 60.0;
  omega= 0.8;
  mx= MX0-1;
  my= MY0-1;
  mz= MZ0-1;
  ndx= NDX0;
  ndy= NDY0;
  ndz= NDZ0;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &npe);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  initcomm(ndx,ndy,ndz);
  it= initmax(mx,my,mz);

  /*
   *    Initializing matrixes
   */
  initmt(mx,it);

  if(id==0){
    printf("Sequential version array size\n");
    printf(" mimax = %d mjmax = %d mkmax = %d\n",MX0,MY0,MZ0);
    printf("Parallel version array size\n");
    printf(" mimax = %d mjmax = %d mkmax = %d\n",MIMAX,MJMAX,MKMAX);
    printf("imax = %d jmax = %d kmax =%d\n",imax,jmax,kmax);
    printf("I-decomp = %d J-decomp = %d K-decomp =%d\n",ndx,ndy,ndz);
  }

  nn= 3;
  if(id==0){
    printf(" Start rehearsal measurement process.\n");
    printf(" Measure the performance in %d times.\n\n",nn);
  }

  MPI_Barrier(MPI_COMM_WORLD);
  cpu0= MPI_Wtime();
  gosa= jacobi(nn);
  cpu1= MPI_Wtime() - cpu0;

  MPI_Allreduce(&cpu1,
                &cpu,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD);

  flop= fflop(mz,my,mx);
  if(id == 0){
    printf(" MFLOPS: %f time(s): %f %e\n\n",
           mflops(nn,cpu,flop),cpu,gosa);
  }
  nn= (int)(target/(cpu/3.0));

  if(id == 0){
    printf(" Now, start the actual measurement process.\n");
    printf(" The loop will be excuted in %d times\n",nn);
    printf(" This will take about one minute.\n");
    printf(" Wait for a while\n\n");
  }

  /*
   *    Start measuring
   */
  MPI_Barrier(MPI_COMM_WORLD);
  cpu0 = MPI_Wtime();
  gosa = jacobi(nn);
  cpu1 = MPI_Wtime() - cpu0;

  MPI_Allreduce(&cpu1,
                &cpu,
                1,
                MPI_DOUBLE,
                MPI_MAX,
                MPI_COMM_WORLD);

  if(id == 0){
    printf("cpu : %f sec.\n", cpu);
    printf("Loop executed for %d times\n",nn);
    printf("Gosa : %e \n",gosa);
    printf("MFLOPS measured : %f\n",mflops(nn,cpu,flop));
    printf("Score based on Pentium III 600MHz : %f\n",
           mflops(nn,cpu,flop)/82.84);
    //todo
    printf("Compression rate: %f \n", 1/(cr/cr_num));
  }

  MPI_Finalize();
  
  return (0);
}

double
fflop(int mx,int my, int mz)
{
  return((double)(mz-2)*(double)(my-2)*(double)(mx-2)*34.0);
}

double
mflops(int nn,double cpu,double flop)
{
  return(flop/cpu*1.e-6*(double)nn);
}

void
initmt(int mx,int it)
{
  int i,j,k;

  for(i=0 ; i<MIMAX ; ++i)
    for(j=0 ; j<MJMAX ; ++j)
      for(k=0 ; k<MKMAX ; ++k){
        a[0][i][j][k]=0.0;
        a[1][i][j][k]=0.0;
        a[2][i][j][k]=0.0;
        a[3][i][j][k]=0.0;
        b[0][i][j][k]=0.0;
        b[1][i][j][k]=0.0;
        b[2][i][j][k]=0.0;
        c[0][i][j][k]=0.0;
        c[1][i][j][k]=0.0;
        c[2][i][j][k]=0.0;
        p[i][j][k]=0.0;
        wrk1[i][j][k]=0.0;
        wrk2[i][j][k]=0.0;
        bnd[i][j][k]=0.0;
      }

  for(i=0 ; i<imax ; ++i)
    for(j=0 ; j<jmax ; ++j)
      for(k=0 ; k<kmax ; ++k){
        a[0][i][j][k]=1.0;
        a[1][i][j][k]=1.0;
        a[2][i][j][k]=1.0;
        a[3][i][j][k]=1.0/6.0;
        b[0][i][j][k]=0.0;
        b[1][i][j][k]=0.0;
        b[2][i][j][k]=0.0;
        c[0][i][j][k]=1.0;
        c[1][i][j][k]=1.0;
        c[2][i][j][k]=1.0;
        p[i][j][k]=(float)((i+it)*(i+it))/(float)((mx-1)*(mx-1));
        wrk1[i][j][k]=0.0;
        wrk2[i][j][k]=0.0;
        bnd[i][j][k]=1.0;
      }
}

float
jacobi(int nn)
{
  int i,j,k,n;
  float gosa,wgosa,s0,ss;

  for(n=0 ; n<nn ; ++n){
    gosa = 0.0;
    wgosa= 0.0;

    for(i=1 ; i<imax-1 ; ++i)
      for(j=1 ; j<jmax-1 ; ++j)
        for(k=1 ; k<kmax-1 ; ++k){
          s0 = a[0][i][j][k] * p[i+1][j  ][k  ]
             + a[1][i][j][k] * p[i  ][j+1][k  ]
             + a[2][i][j][k] * p[i  ][j  ][k+1]
             + b[0][i][j][k] * ( p[i+1][j+1][k  ] - p[i+1][j-1][k  ]
                               - p[i-1][j+1][k  ] + p[i-1][j-1][k  ] )
             + b[1][i][j][k] * ( p[i  ][j+1][k+1] - p[i  ][j-1][k+1]
                               - p[i  ][j+1][k-1] + p[i  ][j-1][k-1] )
             + b[2][i][j][k] * ( p[i+1][j  ][k+1] - p[i-1][j  ][k+1]
                               - p[i+1][j  ][k-1] + p[i-1][j  ][k-1] )
             + c[0][i][j][k] * p[i-1][j  ][k  ]
             + c[1][i][j][k] * p[i  ][j-1][k  ]
             + c[2][i][j][k] * p[i  ][j  ][k-1]
             + wrk1[i][j][k];

          ss = ( s0 * a[3][i][j][k] - p[i][j][k] ) * bnd[i][j][k];
          wgosa += ss*ss;

          wrk2[i][j][k] = p[i][j][k] + omega * ss;
        }

    for(i=1 ; i<imax-1 ; ++i)
      for(j=1 ; j<jmax-1 ; ++j)
        for(k=1 ; k<kmax-1 ; ++k)
          p[i][j][k] = wrk2[i][j][k];

    double start_time = MPI_Wtime();
    sendp(ndx,ndy,ndz);
    double end_time = MPI_Wtime();
    //printf("execution time: %f\n", end_time-start_time);

    MPI_Allreduce(&wgosa,
                  &gosa,
                  1,
                  MPI_FLOAT,
                  MPI_SUM,
                  MPI_COMM_WORLD);
  } /* end n loop */

  return(gosa);
}


void
initcomm(int ndx,int ndy,int ndz)
{
  int  i,j,k,tmp;
  int  ipd[3],idm[3],ir;
  MPI_Comm  icomm;

  if(ndx*ndy*ndz != npe){
    if(id==0){
      printf("Invalid number of PE\n");
      printf("Please check partitioning pattern or number of PE\n");
    }
    MPI_Finalize();
    exit(0);
  }

  icomm= MPI_COMM_WORLD;

  idm[0]= ndx;
  idm[1]= ndy;
  idm[2]= ndz;

  ipd[0]= 0;
  ipd[1]= 0;
  ipd[2]= 0;
  ir= 0;


  MPI_Cart_create(icomm,
                  ndims,
                  idm,
                  ipd,
                  ir,
                  &mpi_comm_cart);
  MPI_Cart_get(mpi_comm_cart,
               ndims,
               idm,
               ipd,
               iop);

  if(ndz > 1){
    MPI_Cart_shift(mpi_comm_cart,
                   2,
                   1,
                   &npz[0],
                   &npz[1]);
  }                     
  if(ndy > 1){
    MPI_Cart_shift(mpi_comm_cart,
                   1,
                   1,
                   &npy[0],
                   &npy[1]);
  }                     
  if(ndx > 1){
    MPI_Cart_shift(mpi_comm_cart,
                   0,
                   1,
                   &npx[0],
                   &npx[1]);
  }                     

}

int
initmax(int mx,int my,int mz)
{
  int  i,tmp,it;
  int  mx1[NDX0+1],my1[NDY0+1],mz1[NDZ0+1];
  int  mx2[NDX0+1],my2[NDY0+1],mz2[NDZ0+1];

  tmp= mx/ndx;
  mx1[0]= 0;
  for(i=1;i<=ndx;i++){
    if(i <= mx%ndx)
      mx1[i]= mx1[i-1] + tmp + 1;
    else
      mx1[i]= mx1[i-1] + tmp;
  }
  tmp= my/ndy;
  my1[0]= 0;
  for(i=1;i<=ndy;i++){
    if(i <= my%ndy)
      my1[i]= my1[i-1] + tmp + 1;
    else
      my1[i]= my1[i-1] + tmp;
  }
  tmp= mz/ndz;
  mz1[0]= 0;
  for(i=1;i<=ndz;i++){
    if(i <= mz%ndz)
      mz1[i]= mz1[i-1] + tmp + 1;
    else
      mz1[i]= mz1[i-1] + tmp;
  }

  for(i=0 ; i<ndx ; i++){
    mx2[i] = mx1[i+1] - mx1[i];
    if(i != 0)     mx2[i] = mx2[i] + 1;
    if(i != ndx-1) mx2[i] = mx2[i] + 1;
  }
  for(i=0 ; i<ndy ; i++){
    my2[i] = my1[i+1] - my1[i];
    if(i != 0)     my2[i] = my2[i] + 1;
    if(i != ndy-1) my2[i] = my2[i] + 1;
  }
  for(i=0 ; i<ndz ; i++){
    mz2[i] = mz1[i+1] - mz1[i];
    if(i != 0)     mz2[i] = mz2[i] + 1;
    if(i != ndz-1) mz2[i] = mz2[i] + 1;
  }

  imax = mx2[iop[0]];
  jmax = my2[iop[1]];
  kmax = mz2[iop[2]];

  if(iop[0] == 0)
    it= mx1[iop[0]];
  else
    it= mx1[iop[0]] - 1;

  if(ndx > 1){
    MPI_Type_vector(jmax,
                    kmax,
                    MKMAX,
                    MPI_FLOAT,
                    &jkvec);
    MPI_Type_commit(&jkvec);
  }                    
  if(ndy > 1){
    MPI_Type_vector(imax,
                    kmax,
                    MJMAX*MKMAX,
                    MPI_FLOAT,
                    &ikvec);
    MPI_Type_commit(&ikvec);
  }                    
  if(ndz > 1){
    MPI_Type_vector(imax*jmax,
                    1,
                    MKMAX,
                    MPI_FLOAT,
                    &ijvec);
    MPI_Type_commit(&ijvec);
  }                    

  return(it);
}

void
sendp(int ndx,int ndy,int ndz)
{
  if(ndz > 1)
    sendp3();

  if(ndy > 1)
    sendp2();

  if(ndx > 1)
    sendp1();
}

void
sendp3()
{
  MPI_Status   st[4];
  MPI_Request  req[4];

  //todo
  MPI_Status   st_plus[12];
  MPI_Request  req_plus[12];

  MPI_Status   st_bitwise[8];
  MPI_Request  req_bitwise[8];  

  if(CT == 5)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npz[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npz[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 3, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 3, kmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*jmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*jmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise(data_small[0], imax*jmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise(data_small[1], imax*jmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npz[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npz[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npz %d %d \n", npz[0], npz[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(imax*jmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*jmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npz[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npz[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npz[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npz[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npz[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npz[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npz[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npz[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise(data_bits_recv[0], data_bytes_recv[0], imax*jmax);
    decompressed_data[1] = myDecompress_bitwise(data_bits_recv[1], data_bytes_recv[1], imax*jmax);
    int pointer = 0;
	int a,b;
    for(a=0; a<imax; a++)
    {
      for(b=0; b<jmax; b++)
      {
        p[a][b][kmax-1] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][b][0] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }

  //revised
  if(CT == 1)
  {
    int array_float_len_send[2] = {0, 0};
    int array_float_len_recv[2] = {0, 0};

    float* array_float_send[2] = {NULL, NULL};
    char* array_char_send[2] = {NULL, NULL};
    int* array_char_displacement_send[2] = {NULL, NULL};

    MPI_Irecv(&array_float_len_recv[0], 1, MPI_INT, npz[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&array_float_len_recv[1], 1, MPI_INT, npz[0], 1, mpi_comm_cart, req+1);
    float* data_0 = transform_3d_array_to_1d_array(p, 3, 1, imax, jmax, kmax);
    float* data_1 = transform_3d_array_to_1d_array(p, 3, kmax-2, imax, jmax, kmax);

    array_float_len_send[0] = myCompress(data_0, &array_float_send[0], &array_char_send[0], &array_char_displacement_send[0], imax*jmax);
    array_float_len_send[1] = myCompress(data_1, &array_float_send[1], &array_char_send[1], &array_char_displacement_send[1], imax*jmax);
    MPI_Isend(&array_float_len_send[0], 1, MPI_INT, npz[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&array_float_len_send[1], 1, MPI_INT, npz[1], 1, mpi_comm_cart, req+3); 
    //printf("send1 %d %d \n", array_float_len_send[0], array_float_len_send[1]);
    MPI_Waitall(4, req, st);

    //printf("recv1 %d %d \n", array_float_len_recv[0], array_float_len_recv[1]);

    float* array_float_recv[2]; 
    char* array_char_recv[2]; 
    int* array_char_displacement_recv[2]; 

    int num_p_recv_0 = array_float_len_recv[0], num_c_recv_0 = imax*jmax - array_float_len_recv[0];
    int num_p_recv_1 = array_float_len_recv[1], num_c_recv_1 = imax*jmax - array_float_len_recv[1];
    array_float_recv[0] = (float*) malloc(sizeof(float)*num_p_recv_0);
    array_char_recv[0] = (char*) malloc(sizeof(char)*num_c_recv_0);
    array_char_displacement_recv[0] = (int*) malloc(sizeof(int)*num_c_recv_0);
    array_float_recv[1] = (float*) malloc(sizeof(float)*num_p_recv_1);
    array_char_recv[1] = (char*) malloc(sizeof(char)*num_c_recv_1);
    array_char_displacement_recv[1] = (int*) malloc(sizeof(int)*num_c_recv_1);
    MPI_Irecv(array_float_recv[0], num_p_recv_0, MPI_FLOAT, npz[1], 2, mpi_comm_cart, req_plus);
    MPI_Irecv(array_char_recv[0], num_c_recv_0, MPI_CHAR, npz[1], 3, mpi_comm_cart, req_plus+1);
    MPI_Irecv(array_char_displacement_recv[0], num_c_recv_0, MPI_INT, npz[1], 4, mpi_comm_cart, req_plus+2);    
    MPI_Irecv(array_float_recv[1], num_p_recv_1, MPI_FLOAT, npz[0], 5, mpi_comm_cart, req_plus+3);
    MPI_Irecv(array_char_recv[1], num_c_recv_1, MPI_CHAR, npz[0], 6, mpi_comm_cart, req_plus+4);
    MPI_Irecv(array_char_displacement_recv[1], num_c_recv_1, MPI_INT, npz[0], 7, mpi_comm_cart, req_plus+5);  
    //mycompress
    int num_p_send_0 = array_float_len_send[0], num_c_send_0 = imax*jmax - array_float_len_send[0];
    int num_p_send_1 = array_float_len_send[1], num_c_send_1 = imax*jmax - array_float_len_send[1];
    MPI_Isend(array_float_send[0], num_p_send_0, MPI_FLOAT, npz[0], 2, mpi_comm_cart, req_plus+6); 
    MPI_Isend(array_char_send[0], num_c_send_0, MPI_CHAR, npz[0], 3, mpi_comm_cart, req_plus+7); 
    MPI_Isend(array_char_displacement_send[0], num_c_send_0, MPI_INT, npz[0], 4, mpi_comm_cart, req_plus+8); 
    MPI_Isend(array_float_send[1], num_p_send_1, MPI_FLOAT, npz[1], 5, mpi_comm_cart, req_plus+9); 
    MPI_Isend(array_char_send[1], num_c_send_1, MPI_CHAR, npz[1], 6, mpi_comm_cart, req_plus+10); 
    MPI_Isend(array_char_displacement_send[1], num_c_send_1, MPI_INT, npz[1], 7, mpi_comm_cart, req_plus+11); 
    MPI_Waitall(12, req_plus, st_plus);

    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + 1.0*((float)num_p_send_0/(imax*jmax));
    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + 1.0*((float)num_p_send_1/(imax*jmax));
    cr_num += 2;

    //calculate bitwise compress ratio
    // printf("%d %d\n", num_p_send_0, num_c_send_0);
    // printf("%d %d\n", num_p_send_1, num_c_send_1);
    // printf("%f \n", calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0));
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax));
    // printf("%f \n", (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax)));
    // printf("%f \n", cr);
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[1], num_p_send_1)*((float)num_p_send_1/(imax*jmax));
    // printf("%f \n", cr);
    // cr_num += 2;
    // printf("%d \n", cr_num);

    float* decompressed_data_0 = myDecompress(array_float_recv[0], array_char_recv[0], array_char_displacement_recv[0], imax*jmax);
    int pointer_0 = 0;
    // for(int a=0; a<imax; a++)
    // {
    //   for(int b=0; b<jmax; b++)
    //   {
    //     p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
    //   }
    // }
    float* decompressed_data_1 = myDecompress(array_float_recv[1], array_char_recv[1], array_char_displacement_recv[1], imax*jmax);
    int pointer_1 = 0;
	int a,b=0;
    for(a=0; a<imax; a++)
    {
      for(b=0; b<jmax; b++)
      {
        p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
        p[a][b][0] = decompressed_data_1[pointer_1++];
      }
    }
  }

  if (CT == 0)
  {
    MPI_Irecv(&p[0][0][kmax-1],
              1,
              ijvec,
              npz[1],
              1,
              mpi_comm_cart,
              req);
    MPI_Irecv(&p[0][0][0],
              1,
              ijvec,
              npz[0],
              2,
              mpi_comm_cart,
              req+1);
  }

  //todo
  // if(CT == 1)
  // {
  //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 3, kmax-2, imax, jmax, kmax);
  // }
  // else if(CT == 2)
  // {
  //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 3, kmax-2, imax, jmax, kmax); 
  // } 
  // else if(CT == 3)    
  // {
  //   cr += calcCompressionRatio_himeno_nolossy_area(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_nolossy_area(p, 3, kmax-2, imax, jmax, kmax);     
  // }
  // else if(CT == 4)    
  // {
  //   cr += calcCompressionRatio_himeno_sz(p, 3, 1, imax, jmax, kmax);
  //   cr += calcCompressionRatio_himeno_sz(p, 3, kmax-2, imax, jmax, kmax);     
  // }
  // cr_num += 2;

  if(CT == 0)
  {
    MPI_Isend(&p[0][0][1],
              1,
              ijvec,
              npz[0],
              1,
              mpi_comm_cart,
              req+2);
    MPI_Isend(&p[0][0][kmax-2],
              1,
              ijvec,
              npz[1],
              2,
              mpi_comm_cart,
              req+3);
    MPI_Waitall(4,
                req,
                st);

    //printf("npz %d %d \n", npz[0], npz[1]);  
  }
}

void
sendp2()
{
  MPI_Status  st[4];
  MPI_Request req[4];

  //todo
  MPI_Status   st_plus[12];
  MPI_Request  req_plus[12];

  MPI_Status   st_bitwise[8];
  MPI_Request  req_bitwise[8];  

  if(CT == 5)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npy[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npy[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 2, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 2, jmax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], imax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], imax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise(data_small[0], imax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise(data_small[1], imax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npy[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npy[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npy %d %d \n", npy[0], npy[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(imax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(imax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npy[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npy[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npy[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npy[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npy[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npy[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npy[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npy[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise(data_bits_recv[0], data_bytes_recv[0], imax*kmax);
    decompressed_data[1] = myDecompress_bitwise(data_bits_recv[1], data_bytes_recv[1], imax*kmax);
    int pointer = 0;
	int a,b;
    for(a=0; a<imax; a++)
    {
      for(b=0; b<kmax; b++)
      {
        p[a][jmax-1][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[a][0][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }  

  if(CT == 1)
  {
    int array_float_len_send[2] = {0, 0};
    int array_float_len_recv[2] = {0, 0};

    float* array_float_send[2] = {NULL, NULL};
    char* array_char_send[2] = {NULL, NULL};
    int* array_char_displacement_send[2] = {NULL, NULL};

    MPI_Irecv(&array_float_len_recv[0], 1, MPI_INT, npy[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&array_float_len_recv[1], 1, MPI_INT, npy[0], 1, mpi_comm_cart, req+1);
    float* data_0 = transform_3d_array_to_1d_array(p, 2, 1, imax, jmax, kmax);
    float* data_1 = transform_3d_array_to_1d_array(p, 2, jmax-2, imax, jmax, kmax);
    array_float_len_send[0] = myCompress(data_0, &array_float_send[0], &array_char_send[0], &array_char_displacement_send[0], imax*kmax);
    array_float_len_send[1] = myCompress(data_1, &array_float_send[1], &array_char_send[1], &array_char_displacement_send[1], imax*kmax);
    MPI_Isend(&array_float_len_send[0], 1, MPI_INT, npy[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&array_float_len_send[1], 1, MPI_INT, npy[1], 1, mpi_comm_cart, req+3); 
    //printf("send1 %d %d \n", array_float_len_send[0], array_float_len_send[1]);
    MPI_Waitall(4, req, st);

    //printf("recv1 %d %d \n", array_float_len_recv[0], array_float_len_recv[1]);

    float* array_float_recv[2]; 
    char* array_char_recv[2]; 
    int* array_char_displacement_recv[2]; 

    int num_p_recv_0 = array_float_len_recv[0], num_c_recv_0 = imax*kmax - array_float_len_recv[0];
    int num_p_recv_1 = array_float_len_recv[1], num_c_recv_1 = imax*kmax - array_float_len_recv[1];
    array_float_recv[0] = (float*) malloc(sizeof(float)*num_p_recv_0);
    array_char_recv[0] = (char*) malloc(sizeof(char)*num_c_recv_0);
    array_char_displacement_recv[0] = (int*) malloc(sizeof(int)*num_c_recv_0);
    array_float_recv[1] = (float*) malloc(sizeof(float)*num_p_recv_1);
    array_char_recv[1] = (char*) malloc(sizeof(char)*num_c_recv_1);
    array_char_displacement_recv[1] = (int*) malloc(sizeof(int)*num_c_recv_1);
    MPI_Irecv(array_float_recv[0], num_p_recv_0, MPI_FLOAT, npy[1], 2, mpi_comm_cart, req_plus);
    MPI_Irecv(array_char_recv[0], num_c_recv_0, MPI_CHAR, npy[1], 3, mpi_comm_cart, req_plus+1);
    MPI_Irecv(array_char_displacement_recv[0], num_c_recv_0, MPI_INT, npy[1], 4, mpi_comm_cart, req_plus+2);    
    MPI_Irecv(array_float_recv[1], num_p_recv_1, MPI_FLOAT, npy[0], 5, mpi_comm_cart, req_plus+3);
    MPI_Irecv(array_char_recv[1], num_c_recv_1, MPI_CHAR, npy[0], 6, mpi_comm_cart, req_plus+4);
    MPI_Irecv(array_char_displacement_recv[1], num_c_recv_1, MPI_INT, npy[0], 7, mpi_comm_cart, req_plus+5);  
    //mycompress
    int num_p_send_0 = array_float_len_send[0], num_c_send_0 = imax*kmax - array_float_len_send[0];
    int num_p_send_1 = array_float_len_send[1], num_c_send_1 = imax*kmax - array_float_len_send[1];
    MPI_Isend(array_float_send[0], num_p_send_0, MPI_FLOAT, npy[0], 2, mpi_comm_cart, req_plus+6); 
    MPI_Isend(array_char_send[0], num_c_send_0, MPI_CHAR, npy[0], 3, mpi_comm_cart, req_plus+7); 
    MPI_Isend(array_char_displacement_send[0], num_c_send_0, MPI_INT, npy[0], 4, mpi_comm_cart, req_plus+8); 
    MPI_Isend(array_float_send[1], num_p_send_1, MPI_FLOAT, npy[1], 5, mpi_comm_cart, req_plus+9); 
    MPI_Isend(array_char_send[1], num_c_send_1, MPI_CHAR, npy[1], 6, mpi_comm_cart, req_plus+10); 
    MPI_Isend(array_char_displacement_send[1], num_c_send_1, MPI_INT, npy[1], 7, mpi_comm_cart, req_plus+11); 
    MPI_Waitall(12, req_plus, st_plus);

    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*kmax)) + 1.0*((float)num_p_send_0/(imax*kmax));
    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*kmax)) + 1.0*((float)num_p_send_1/(imax*kmax));
    cr_num += 2;

    //calculate bitwise compress ratio
    // printf("%d %d\n", num_p_send_0, num_c_send_0);
    // printf("%d %d\n", num_p_send_1, num_c_send_1);
    // printf("%f \n", calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0));
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax));
    // printf("%f \n", (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax)));
    // printf("%f \n", cr);
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[1], num_p_send_1)*((float)num_p_send_1/(imax*jmax));
    // printf("%f \n", cr);
    // cr_num += 2;
    // printf("%d \n", cr_num);

    float* decompressed_data_0 = myDecompress(array_float_recv[0], array_char_recv[0], array_char_displacement_recv[0], imax*kmax);
    int pointer_0 = 0;
    // for(int a=0; a<imax; a++)
    // {
    //   for(int b=0; b<jmax; b++)
    //   {
    //     p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
    //   }
    // }
    float* decompressed_data_1 = myDecompress(array_float_recv[1], array_char_recv[1], array_char_displacement_recv[1], imax*kmax);
    int pointer_1 = 0;
	int a,b;
    for(a=0; a<imax; a++)
    {
      for(b=0; b<kmax; b++)
      {
        p[a][jmax-1][b] = decompressed_data_0[pointer_0++];
        p[a][0][b] = decompressed_data_1[pointer_1++];
      }
    }
  }  

  if(CT == 0)
  {
    MPI_Irecv(&p[0][jmax-1][0],
              1,
              ikvec,
              npy[1],
              1,
              mpi_comm_cart,
              req);
    MPI_Irecv(&p[0][0][0],
              1,
              ikvec,
              npy[0],
              2,
              mpi_comm_cart,
              req+1);
    //todo
    // if (CT == 1)
    // {
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 2, jmax-2, imax, jmax, kmax); 
    // }
    // else if(CT == 2)
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 2, jmax-2, imax, jmax, kmax); 
    // } 
    // else if(CT == 3)    
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 2, jmax-2, imax, jmax, kmax);     
    // }   
    // else if(CT == 4)    
    // {
    //   cr += calcCompressionRatio_himeno_sz(p, 2, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_sz(p, 2, jmax-2, imax, jmax, kmax);     
    // }      
    // cr_num += 2;         
    MPI_Isend(&p[0][1][0],
              1,
              ikvec,
              npy[0],
              1,
              mpi_comm_cart,
              req+2);
    MPI_Isend(&p[0][jmax-2][0],
              1,
              ikvec,
              npy[1],
              2,
              mpi_comm_cart,
              req+3);

    MPI_Waitall(4,
                req,
                st);
  }


}


void
sendp1()
{
  MPI_Status  st[4];
  MPI_Request req[4];

  //todo
  MPI_Status   st_plus[12];
  MPI_Request  req_plus[12];

  MPI_Status   st_bitwise[8];
  MPI_Request  req_bitwise[8];  

  if(CT == 5)
  {
    int data_bytes_send[2] = {0, 0};
    int data_bytes_recv[2] = {0, 0};

    MPI_Irecv(&data_bytes_recv[0], 1, MPI_INT, npx[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&data_bytes_recv[1], 1, MPI_INT, npx[0], 1, mpi_comm_cart, req+1);
    
    float* data[2];
    data[0] = transform_3d_array_to_1d_array(p, 1, 1, imax, jmax, kmax);
    data[1] = transform_3d_array_to_1d_array(p, 1, imax-2, imax, jmax, kmax);
    float* data_small[2] = {NULL, NULL};
    float data_min_send[2];
    data_min_send[0] = toSmallDataset_float(data[0], &data_small[0], jmax*kmax);        
    data_min_send[1] = toSmallDataset_float(data[1], &data_small[1], jmax*kmax);    

    unsigned char* data_bits_send[2] = {NULL, NULL};

    int data_pos[2] = {8, 8}; //position of filled bit in last byte --> 87654321

    myCompress_bitwise(data_small[0], jmax*kmax, &data_bits_send[0], &data_bytes_send[0], &data_pos[0]);
    myCompress_bitwise(data_small[1], jmax*kmax, &data_bits_send[1], &data_bytes_send[1], &data_pos[1]);

    MPI_Isend(&data_bytes_send[0], 1, MPI_INT, npx[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&data_bytes_send[1], 1, MPI_INT, npx[1], 1, mpi_comm_cart, req+3); 
    MPI_Waitall(4, req, st);  

    // printf("npx %d %d \n", npx[0], npx[1]);  
    // printf("send %d %d \n", data_bytes_send[0], data_bytes_send[1]);
    // printf("recv %d %d \n", data_bytes_recv[0], data_bytes_recv[1]);

    cr += data_bytes_send[0]*8.0/(jmax*kmax*sizeof(float)*8);
    cr += data_bytes_send[1]*8.0/(jmax*kmax*sizeof(float)*8);
    cr_num += 2;

    unsigned char* data_bits_recv[2];
    float data_min_recv[2];

    data_bits_recv[0] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[0]);
    data_bits_recv[1] = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_recv[1]);
    MPI_Irecv(&data_min_recv[0], 1, MPI_FLOAT, npx[1], 2, mpi_comm_cart, req_bitwise);
    MPI_Irecv(data_bits_recv[0], data_bytes_recv[0], MPI_UNSIGNED_CHAR, npx[1], 3, mpi_comm_cart, req_bitwise+1);
    MPI_Irecv(&data_min_recv[1], 1, MPI_FLOAT, npx[0], 4, mpi_comm_cart, req_bitwise+2);
    MPI_Irecv(data_bits_recv[1], data_bytes_recv[1], MPI_UNSIGNED_CHAR, npx[0], 5, mpi_comm_cart, req_bitwise+3);  
    //bitwise compress
    MPI_Isend(&data_min_send[0], 1, MPI_FLOAT, npx[0], 2, mpi_comm_cart, req_bitwise+4); 
    MPI_Isend(data_bits_send[0], data_bytes_send[0], MPI_UNSIGNED_CHAR, npx[0], 3, mpi_comm_cart, req_bitwise+5); 
    MPI_Isend(&data_min_send[1], 1, MPI_FLOAT, npx[1], 4, mpi_comm_cart, req_bitwise+6); 
    MPI_Isend(data_bits_send[1], data_bytes_send[1], MPI_UNSIGNED_CHAR, npx[1], 5, mpi_comm_cart, req_bitwise+7); 
    MPI_Waitall(8, req_bitwise, st_bitwise);    

    float* decompressed_data[2];
    decompressed_data[0] = myDecompress_bitwise(data_bits_recv[0], data_bytes_recv[0], jmax*kmax);
    decompressed_data[1] = myDecompress_bitwise(data_bits_recv[1], data_bytes_recv[1], jmax*kmax);
    int pointer = 0;
	int a,b;
    for(a=0; a<jmax; a++)
    {
      for(b=0; b<kmax; b++)
      {
        p[imax-1][a][b] = decompressed_data[0][pointer] + data_min_recv[0];
        p[0][a][b] = decompressed_data[1][pointer] + data_min_recv[1];
        pointer++;
      }
    }
  }  

  if(CT == 1)
  {
    int array_float_len_send[2] = {0, 0};
    int array_float_len_recv[2] = {0, 0};

    float* array_float_send[2] = {NULL, NULL};
    char* array_char_send[2] = {NULL, NULL};
    int* array_char_displacement_send[2] = {NULL, NULL};

    MPI_Irecv(&array_float_len_recv[0], 1, MPI_INT, npx[1], 0, mpi_comm_cart, req);
    MPI_Irecv(&array_float_len_recv[1], 1, MPI_INT, npx[0], 1, mpi_comm_cart, req+1);
    float* data_0 = transform_3d_array_to_1d_array(p, 1, 1, imax, jmax, kmax);
    float* data_1 = transform_3d_array_to_1d_array(p, 1, imax-2, imax, jmax, kmax);
    array_float_len_send[0] = myCompress(data_0, &array_float_send[0], &array_char_send[0], &array_char_displacement_send[0], jmax*kmax);
    array_float_len_send[1] = myCompress(data_1, &array_float_send[1], &array_char_send[1], &array_char_displacement_send[1], jmax*kmax);
    MPI_Isend(&array_float_len_send[0], 1, MPI_INT, npx[0], 0, mpi_comm_cart, req+2); 
    MPI_Isend(&array_float_len_send[1], 1, MPI_INT, npx[1], 1, mpi_comm_cart, req+3); 
    //printf("send1 %d %d \n", array_float_len_send[0], array_float_len_send[1]);
    MPI_Waitall(4, req, st);

    //printf("recv1 %d %d \n", array_float_len_recv[0], array_float_len_recv[1]);

    float* array_float_recv[2]; 
    char* array_char_recv[2]; 
    int* array_char_displacement_recv[2]; 

    int num_p_recv_0 = array_float_len_recv[0], num_c_recv_0 = jmax*kmax - array_float_len_recv[0];
    int num_p_recv_1 = array_float_len_recv[1], num_c_recv_1 = jmax*kmax - array_float_len_recv[1];
    array_float_recv[0] = (float*) malloc(sizeof(float)*num_p_recv_0);
    array_char_recv[0] = (char*) malloc(sizeof(char)*num_c_recv_0);
    array_char_displacement_recv[0] = (int*) malloc(sizeof(int)*num_c_recv_0);
    array_float_recv[1] = (float*) malloc(sizeof(float)*num_p_recv_1);
    array_char_recv[1] = (char*) malloc(sizeof(char)*num_c_recv_1);
    array_char_displacement_recv[1] = (int*) malloc(sizeof(int)*num_c_recv_1);
    MPI_Irecv(array_float_recv[0], num_p_recv_0, MPI_FLOAT, npx[1], 2, mpi_comm_cart, req_plus);
    MPI_Irecv(array_char_recv[0], num_c_recv_0, MPI_CHAR, npx[1], 3, mpi_comm_cart, req_plus+1);
    MPI_Irecv(array_char_displacement_recv[0], num_c_recv_0, MPI_INT, npx[1], 4, mpi_comm_cart, req_plus+2);    
    MPI_Irecv(array_float_recv[1], num_p_recv_1, MPI_FLOAT, npx[0], 5, mpi_comm_cart, req_plus+3);
    MPI_Irecv(array_char_recv[1], num_c_recv_1, MPI_CHAR, npx[0], 6, mpi_comm_cart, req_plus+4);
    MPI_Irecv(array_char_displacement_recv[1], num_c_recv_1, MPI_INT, npx[0], 7, mpi_comm_cart, req_plus+5);  
    //mycompress
    int num_p_send_0 = array_float_len_send[0], num_c_send_0 = jmax*kmax - array_float_len_send[0];
    int num_p_send_1 = array_float_len_send[1], num_c_send_1 = jmax*kmax - array_float_len_send[1];
    MPI_Isend(array_float_send[0], num_p_send_0, MPI_FLOAT, npx[0], 2, mpi_comm_cart, req_plus+6); 
    MPI_Isend(array_char_send[0], num_c_send_0, MPI_CHAR, npx[0], 3, mpi_comm_cart, req_plus+7); 
    MPI_Isend(array_char_displacement_send[0], num_c_send_0, MPI_INT, npx[0], 4, mpi_comm_cart, req_plus+8); 
    MPI_Isend(array_float_send[1], num_p_send_1, MPI_FLOAT, npx[1], 5, mpi_comm_cart, req_plus+9); 
    MPI_Isend(array_char_send[1], num_c_send_1, MPI_CHAR, npx[1], 6, mpi_comm_cart, req_plus+10); 
    MPI_Isend(array_char_displacement_send[1], num_c_send_1, MPI_INT, npx[1], 7, mpi_comm_cart, req_plus+11); 
    MPI_Waitall(12, req_plus, st_plus);

    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_0/(jmax*kmax)) + 1.0*((float)num_p_send_0/(jmax*kmax));
    cr += (8.0/(sizeof(float)*8))*((float)num_c_send_1/(jmax*kmax)) + 1.0*((float)num_p_send_1/(jmax*kmax));
    cr_num += 2;

    //calculate bitwise compress ratio
    // printf("%d %d\n", num_p_send_0, num_c_send_0);
    // printf("%d %d\n", num_p_send_1, num_c_send_1);
    // printf("%f \n", calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0));
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax));
    // printf("%f \n", (3.0/(sizeof(float)*8))*((float)num_c_send_0/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[0], num_p_send_0)*((float)num_p_send_0/(imax*jmax)));
    // printf("%f \n", cr);
    // cr += (3.0/(sizeof(float)*8))*((float)num_c_send_1/(imax*jmax)) + calCompressRatio_bitwise_float(array_float_send[1], num_p_send_1)*((float)num_p_send_1/(imax*jmax));
    // printf("%f \n", cr);
    // cr_num += 2;
    // printf("%d \n", cr_num);

    float* decompressed_data_0 = myDecompress(array_float_recv[0], array_char_recv[0], array_char_displacement_recv[0], jmax*kmax);
    int pointer_0 = 0;
    // for(int a=0; a<imax; a++)
    // {
    //   for(int b=0; b<jmax; b++)
    //   {
    //     p[a][b][kmax-1] = decompressed_data_0[pointer_0++];
    //   }
    // }
    float* decompressed_data_1 = myDecompress(array_float_recv[1], array_char_recv[1], array_char_displacement_recv[1], jmax*kmax);
    int pointer_1 = 0;
	int a,b;
    for(a=0; a<jmax; a++)
    {
      for(b=0; b<kmax; b++)
      {
        p[imax-1][a][b] = decompressed_data_0[pointer_0++];
        p[0][a][b] = decompressed_data_1[pointer_1++];
      }
    }
  } 

  if(CT == 0)
  {
    MPI_Irecv(&p[imax-1][0][0],
              1,
              jkvec,
              npx[1],
              1,
              mpi_comm_cart,
              req);
    MPI_Irecv(&p[0][0][0],
              1,
              jkvec,
              npx[0],
              2,
              mpi_comm_cart,
              req+1);
    //todo
    // if (CT == 1)
    // {
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_ij_ik_jk(p, 1, imax-2, imax, jmax, kmax); 
    // }
    // else if(CT == 2)
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_performance(p, 1, imax-2, imax, jmax, kmax); 
    // } 
    // else if(CT == 3)    
    // {
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_nolossy_area(p, 1, imax-2, imax, jmax, kmax);     
    // }    
    // else if(CT == 4)    
    // {
    //   cr += calcCompressionRatio_himeno_sz(p, 1, 1, imax, jmax, kmax);
    //   cr += calcCompressionRatio_himeno_sz(p, 1, imax-2, imax, jmax, kmax);     
    // }   
    // cr_num += 2;               
    MPI_Isend(&p[1][0][0],
              1,
              jkvec,
              npx[0],
              1,
              mpi_comm_cart,
              req+2);
    MPI_Isend(&p[imax-2][0][0],
              1,
              jkvec,
              npx[1],
              2,
              mpi_comm_cart,
              req+3);

    MPI_Waitall(4,
                req,
                st);
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

  int i;
  for(i=0; i<bytes; i++)
  {
    //if(i == bytes - 1 && pos != 8) min_shift = pos;

    int j;
    for (j=7; j>=min_shift; j--) //each bit of byte
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
	  int i;
          for(i=1; i<12; i++)
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
		int i;
        for(i=bits_num+1; i< sizeof(double)*8; i++)
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
  int i;
  for(i=0; i<bytes; i++)
  {
    //if(i == bytes - 1 && pos != 8) min_shift = pos;
    int j;
    for (j=7; j>=min_shift; j--) //each bit of byte
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
	  int i;
          for(i=1; i<9; i++)
          {
            expo_value += (bits[i]-'0')*pow(2,8-i);
          }
          expo_value -= 127;

          int mantissa_bits_within_error_bound = absErrBound_binary + expo_value;
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
	int i;
        for(i=bits_num+1; i< sizeof(float)*8; i++)
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

	int n;
  for(n=0; n<num; n++)
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
  int n;
  for(n=0; n<num; n++)
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
        compress_bitwise_float(real_value, data_bits, bytes, pos);            
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
  int i;
  for(i=1; i<12; i++)
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
  //int i;
  for(i=0; i<bits_after_compress; i++)
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
  int i;
  for(i=1; i<9; i++)
  {
    expo_value += (float_arr[i]-'0')*pow(2,8-i);
  }
  expo_value -= 127;

  int mantissa_bits_within_error_bound = absErrBound_binary + expo_value;

  if(mantissa_bits_within_error_bound > 23) //23 mantissa part of float (52 in the case of double)
  {
    mantissa_bits_within_error_bound = 23;
  }
  else if(mantissa_bits_within_error_bound < 0)
  {
    mantissa_bits_within_error_bound = 0;
  }

  int bits_after_compress = 1+8+mantissa_bits_within_error_bound;  
  //int i;
  for(i=0; i<bits_after_compress; i++)
  {
    add_bit_to_bytes(data_bits, bytes, pos, float_arr[i]-'0');
  }
}

double toSmallDataset_double(double data[], double** data_small, int num)
{
  *data_small = malloc(sizeof(double) * num);
  double min = data[0];
  int i;
  for(i=1; i<num; i++)
  {
    if(data[i]<min)
    {
      min = data[i];
    }
  }
  //int i;
  for(i=0; i<num; i++)
  {
    (*data_small)[i] = data[i] - min;
  }

  return min;
}

float toSmallDataset_float(float data[], float** data_small, int num)
{
  *data_small = malloc(sizeof(float) * num);
  float min = data[0];
  int i;
  for(i=1; i<num; i++)
  {
    if(data[i]<min)
    {
      min = data[i];
    }
  }
  //  int i;
  for(i=0; i<num; i++)
  {
    (*data_small)[i] = data[i] - min;
  }

  return min;
}

float calCompressRatio_bitwise_double2(float data[], int num)
{
  int bits_after_compress = 0;
  int n;
  for(n=0; n<num; n++)
  {
    double double10 = data[n];
    char double_arr[64+1];
    doubletostr(&double10, double_arr);
    //printf("%s \n", double_arr); 

    int expo_value = 0;
    //printf("%f \n", double10); 
    int i;
    for(i=1; i<12; i++)
    {
      expo_value += (double_arr[i]-'0')*pow(2,11-i);
      //printf("%d ", double_arr[i]-'0'); 
    }
    expo_value -= 1023;

    //printf("%d ", expo_value); 

    int mantissa_bits_within_error_bound = absErrBound_binary + expo_value;

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
  int n;
  for(n=0; n<num; n++)
  {
    double double10 = data[n];
    char double_arr[64+1];
    doubletostr(&double10, double_arr);
    //printf("%s \n", double_arr); 

    int expo_value = 0;
    //printf("%f \n", double10); 
    int i;
    for(i=1; i<12; i++)
    {
      expo_value += (double_arr[i]-'0')*pow(2,11-i);
      //printf("%d ", double_arr[i]-'0'); 
    }
    expo_value -= 1023;

    //printf("%d ", expo_value); 

    int mantissa_bits_within_error_bound = absErrBound_binary + expo_value;

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
  int n;
  for(n=0; n<num; n++)
  {
    float float10 = data[n];
    char float_arr[32+1];
    floattostr(&float10, float_arr);
    //printf("%s \n", float_arr); 

    int expo_value = 0;
    int i;
    for(i=1; i<9; i++)
    {
      expo_value += (float_arr[i]-'0')*pow(2,8-i);
      //printf("%d ", float_arr[i]-'0'); 
    }
    expo_value -= 127;

    //printf("%d ", expo_value); 

    int mantissa_bits_within_error_bound = absErrBound_binary + expo_value;

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
  int a,b;
  for(a=0; a<A; a++)
  {
    for(b=0; b<B; b++)
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
  int i;
  for(i=0; i<num; i++)
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

  int n;
  for(n=0; n<num; n++)
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
  int i;
  for(i=0; i<num; i++)
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

  int n;
  for(n=0; n<num; n++)
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
  int a,b;
  for(a=0; a<A; a++)
  {
    for(b=0; b<B; b++)
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
  int a,b;
  for(a=0; a<A; a++)
  {
    for(b=0; b<B; b++)
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

	  int i;
          for(i=1;i<9;i++) //1-9 exponential part of float (1-12 in the case of double)
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
  int a,b;
  for(a=0; a<A; a++)
  {
    for(b=0; b<B; b++)
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
	int i;
        for(i=1;i<sizeof(float)*8;i++)
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
  int a,b;
  for(a=0; a<A; a++)
  {
    for(b=0; b<B; b++)
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
	int i;
        for(i=1;i<sizeof(float)*8;i++)
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

  int n;
  for(n=0; n<num; n++)
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
	int i;
        for(i=1;i<9;i++) //1-9 exponential part of float (1-12 in the case of double)
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

  int n;
  for(n=0; n<num; n++)
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
      int i;
      for(i=1;i<sizeof(float)*8;i++)
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

	int n;
  for(n=0; n<num; n++)
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
	int i;
      for(i=1;i<sizeof(float)*8;i++)
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

	int n;
  for(n=0; n<num; n++)
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

		int i;
        for(i=1;i<12;i++) //1-9 exponential part of float (1-12 in the case of double)
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

	int n;
  for(n=0; n<num; n++)
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
		int i;
      for(i=1;i<sizeof(double)*8;i++)
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

  int n;
  for(n=0; n<num; n++)
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
	int i;
      for(i=1;i<sizeof(double)*8;i++)
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
	int i;
    for(i=0;i<32;i++)
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
	int i=9;
    for(i=0;i<64;i++)
    {
        bin[i] = (*f)&(t<<63-i)?1:0;
    }
}

//100.0 --> 01000010110010000000000000000000
//str should have at least 33 byte.
void floattostr(float* a, char* str){
	unsigned int c;
	c= ((unsigned int*)a)[0]; 
	int i;
	for(i=0;i<32;i++){
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
	int i;
	for(i=0;i<64;i++){
		str[63-i]=(char)(c&1)+'0';
		c>>=1;
	}
	str[64] = '\0';
}

//01000010110010000000000000000000 --> 100.0
float strtofloat(char * str){
	unsigned int flt = 0;
	int i;
	for(i=0;i<31;i++){
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
	int i=0;
	for(i=0;i<63;i++){
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
    float *arr = malloc(sizeof(float) * count);
 
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
    double *arr = malloc(sizeof(double) * count);
 
    fread(arr, sizeof(double), count, fp);
    fclose(fp);

    return arr;
}

void readfrombinary_writetotxt_float(const char *binaryfile, const char *txtfile, int count)
{
    FILE *fp;
    fp = fopen(binaryfile, "rb");
    assert(fp);

    float *arr = malloc(sizeof(float) * count);
    fread(arr, sizeof(float), count, fp);
    fclose(fp);

    fp = fopen(txtfile, "w");
    assert(fp);
     
	int i;
    for(i=0; i<count; i++)
    {
      fprintf(fp, "%f\n", arr[i]);
    }
    fclose(fp);
}

void readfrombinary_writetotxt_double(const char *binaryfile, const char *txtfile, int count)
{
    FILE *fp;
    fp = fopen(binaryfile, "rb");
    assert(fp);

    double *arr = malloc(sizeof(double) * count);
    fread(arr, sizeof(double), count, fp);
    fclose(fp);

    fp = fopen(txtfile, "w");
    assert(fp);

    int i;
    for(i=0; i<count; i++)
    {
      fprintf(fp, "%lf\n", arr[i]);
    }
    fclose(fp);
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
