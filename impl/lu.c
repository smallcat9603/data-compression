#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <math.h>
#include <unistd.h>
#include "param.h"
#include "dataCompression.h"

#define ln() putchar('\n')
// #define GENERIC_TAG (0)

//todo
struct vector
{
  double* p_data; //precise data
  char* c_data; //compressed data
  int* disp; //displacement of compressed data
};

double *gen_mx (size_t dim);
double *gen_row(size_t dim);
double *gen_row_ref (size_t dim, size_t ref);
void print_mx (double *M, size_t dim, size_t sep);
void forw_elim(double **origin, double *master_row, size_t dim);
void U_print (double *M, int dim);
void L_print (double *M, int dim);

int main(int argc, char *argv[])
{
   srand(time(NULL));
   const int root_p = 0;
   int mx_size = 0, p, id;
   if (argc < 2) {
      printf("Matrix size missing in the arguments\n");
      return EXIT_FAILURE;
   }
   mx_size = atol(argv[1]);
   double *A = gen_mx(mx_size);

   MPI_Init(NULL, NULL);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   MPI_Comm_rank(MPI_COMM_WORLD, &id);

   if (id == root_p) {
      printf("[A]\n");
      print_mx(A, mx_size * mx_size, mx_size);
   }   

   int i, j, tmp_size = mx_size - 1, diag_ref = 0;
   double start = MPI_Wtime();

   double gosa = 0;
   float compress_ratio = 0;
   float sz_comp_ratio = 0;
   float nolossy_performance = 0;
   float nolossy_area = 0;  
   uint32_t crc = 0;
   uint32_t crc_check = 0;
   unsigned char crc_ok = 'y';
   unsigned char* crc_ok_recv = NULL;   
   int resent = 0;
   srand((unsigned)time(NULL));     
   for (i = 0; i < tmp_size; i++, diag_ref++) {
      double *diag_row = &A[diag_ref * mx_size + diag_ref];
      for (j = diag_ref + 1; j < mx_size; j++) {
         if (j % p == id) {
            double *save = &A[j * mx_size + diag_ref];
            //printf("[%d] ", id);
            //print_mx(save, mx_size - diag_ref, mx_size - diag_ref);
            forw_elim(&save, diag_row, mx_size - diag_ref);
         }
      }

      for (j = diag_ref + 1; j < mx_size; j++) {
         double *save = &A[j * mx_size + diag_ref];

         int size = mx_size - diag_ref;   
         int root = j % p;  
         if(CT == 8)
         {
            MPI_Bcast_bitwise_crc(save, size, root, id, p, &compress_ratio, &gosa, &resent);

            // int data_bytes = 0;
            // double min = 0;

            // unsigned char* data_bits = NULL;
            
            // if(id == root)
            // {
            //    // sz_comp_ratio += calcCompressionRatio_sz_double(save, size);
            //    // nolossy_performance += calcCompressionRatio_nolossy_performance_double(save, size);
            //    // nolossy_area += calcCompressionRatio_nolossy_area_double(save, size);

            //    //mycommpress
            //    double* small = NULL;
            //    min = toSmallDataset_double(save, &small, size);

            //    int data_pos = 8; //position of filled bit in last byte --> 87654321

            //    myCompress_bitwise_double(small, size, &data_bits, &data_bytes, &data_pos);	
            //    crc = do_crc32(data_bits, data_bytes);		
            // }

            // MPI_Bcast(&data_bytes, 1, MPI_INT, root, MPI_COMM_WORLD);
            // MPI_Bcast(&min, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
            // compress_ratio += data_bytes*8.0/(size*sizeof(double)*8);
         
            // if(id != root)
            // {
            //       data_bits = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes);
            // }
            // MPI_Bcast(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);
            // MPI_Bcast(&crc, 1, MPI_UNSIGNED, root, MPI_COMM_WORLD);

            // if(id != root)
            // {
            //       crc_check = do_crc32(data_bits, data_bytes);

            //       if(BER > 0)
            //       {
            //          double ber = BER;
            //          uint64_t to = 1/ber;
            //          uint64_t r = get_random_int(0, to);
            //          if(r < data_bytes * 8)
            //          {
            //             crc_check = 0;
            //          }
            //       }
                  
            //       if (crc == crc_check)
            //       {
            //          // printf("CRC passed\n");
            //          crc_ok = 'y';
            //       }  
            //       else
            //       {
            //          // printf("CRC NOT passed\n");
            //          crc_ok = 'n';
            //       }            
            // }
            // else
            // {
            //       crc_ok_recv = (unsigned char *)malloc(p*1*sizeof(unsigned char));
            // }

            // MPI_Gather(&crc_ok, 1, MPI_UNSIGNED_CHAR, crc_ok_recv, 1, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);
            
            // if(id == root)
            // {
            //       for(int i = 0; i < p; i++)
            //       {
            //          if(i != root && crc_ok_recv[i] == 'n')
            //          {
            //             MPI_Send(data_bits, data_bytes, MPI_UNSIGNED_CHAR, i, i, MPI_COMM_WORLD);
            //             resent++;
            //          }
            //       }
            // }
            // else if(crc_ok == 'n')
            // {
            //       MPI_Recv(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // }            

            // double* decompressed_data = myDecompress_bitwise_double(data_bits, data_bytes, size);

            // double gs = 0;
            // for(int i=0; i<size; i++)
            // {
            //       if(id == root)
            //       {
            //          gs += fabs(decompressed_data[i] + min - save[i]);
            //       }
            //       else
            //       {
            //          save[i] = decompressed_data[i] + min;
            //       }
            // }
            // gosa += gs/size;

            // //todo
            // free(data_bits);
         }
         else if(CT == 7)
         {
            int data_bytes = 0;
            double min = 0;

            unsigned char* data_bits = NULL;

            int type = 0;
            double medium = 0;
            
            if(id == root)
            {
               //mycommpress
               double* small = NULL;
               min = toSmallDataset_double(save, &small, size);

               int data_pos = 8; //position of filled bit in last byte --> 87654321

               medium = med_dataset_double(small, size, &type);
               char double_arr[64+1];
               doubletostr(&medium, double_arr);
               char mask[1+11+8];
               strncpy(mask, double_arr, 1+11+8);			

               myCompress_bitwise_double_mask(small, size, &data_bits, &data_bytes, &data_pos, type, mask);			
            }

            MPI_Bcast(&data_bytes, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(&min, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
            compress_ratio += data_bytes*8.0/(size*sizeof(double)*8);
         
            if(id != root)
            {
                  data_bits = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes);
            }
            MPI_Bcast(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);

            MPI_Bcast(&medium, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
            MPI_Bcast(&type, 1, MPI_INT, root, MPI_COMM_WORLD);

            char double_arr_recv[64+1];
            doubletostr(&medium, double_arr_recv);
            char mask_recv[1+11+8];
            strncpy(mask_recv, double_arr_recv, 1+11+8);			

            double* decompressed_data = myDecompress_bitwise_double_mask(data_bits, data_bytes, size, type, mask_recv);

            double gs = 0;
            for(int i=0; i<size; i++)
            {
                  if(id == root)
                  {
                     gs += fabs(decompressed_data[i] + min - save[i]);
                  }
                  else
                  {
                     save[i] = decompressed_data[i] + min;
                  }
            }  
            gosa += gs/size;

            //todo
            free(data_bits);
         }
         else if(CT == 6)
         {
            int data_bytes = 0;
            double min = 0;

            unsigned char* data_bits = NULL;
            
            if(id == root)
            {
                  //mycommpress
                  double* small = NULL;
                  min = toSmallDataset_double(save, &small, size);

                  int data_pos = 8; //position of filled bit in last byte --> 87654321

                  myCompress_bitwise_double_np(small, size, &data_bits, &data_bytes, &data_pos);			
            }

            MPI_Bcast(&data_bytes, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(&min, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
            compress_ratio += data_bytes*8.0/(size*sizeof(double)*8);
         
            if(id != root)
            {
                  data_bits = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes);
            }
            MPI_Bcast(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);

            double* decompressed_data = myDecompress_bitwise_double_np(data_bits, data_bytes, size);

            double gs = 0;
            for(int i=0; i<size; i++)
            {
                  if(id == root)
                  {
                     gs += fabs(decompressed_data[i] + min - save[i]);
                  }
                  else
                  {
                     save[i] = decompressed_data[i] + min;
                  }
            }
            gosa += gs/size;        

            //todo
            free(data_bits);
         }
         else if(CT == 5)
         {
            int data_bytes = 0;
            double min = 0;

            unsigned char* data_bits = NULL;
            
            if(id == root)
            {
               // sz_comp_ratio += calcCompressionRatio_sz_double(save, size);
               // nolossy_performance += calcCompressionRatio_nolossy_performance_double(save, size);
               // nolossy_area += calcCompressionRatio_nolossy_area_double(save, size);

               //mycommpress
               double* small = NULL;
               min = toSmallDataset_double(save, &small, size);

               int data_pos = 8; //position of filled bit in last byte --> 87654321

               myCompress_bitwise_double(small, size, &data_bits, &data_bytes, &data_pos);			
            }

            MPI_Bcast(&data_bytes, 1, MPI_INT, root, MPI_COMM_WORLD);
            MPI_Bcast(&min, 1, MPI_DOUBLE, root, MPI_COMM_WORLD);
            compress_ratio += data_bytes*8.0/(size*sizeof(double)*8);
         
            if(id != root)
            {
                  data_bits = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes);
            }
            MPI_Bcast(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);

            double* decompressed_data = myDecompress_bitwise_double(data_bits, data_bytes, size);

            double gs = 0;
            for(int i=0; i<size; i++)
            {
                  if(id == root)
                  {
                     gs += fabs(decompressed_data[i] + min - save[i]);
                  }
                  else
                  {
                     save[i] = decompressed_data[i] + min;
                  }
            }
            gosa += gs/size;

            //todo
            free(data_bits);
         }
         else if(CT == 4)
         {
            int data_bytes = 0;

            unsigned char* data_bits = NULL;
            
            if(id == root)
            {
               char binfile[64];
               sprintf(binfile, "dataset/id%d.dat", root);
               writetobinary_double(binfile, save, size); //.txt --> .dat
               char sz_comp_cmd[64];
               sprintf(sz_comp_cmd, "%s%g%sdataset/id%d%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, root, sz_comp_cmd_suffix2, size);
               //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
               int iret_comp = system(sz_comp_cmd); //.dat --> .dat.sz
               char binfile_sz[64];
               sprintf(binfile_sz, "dataset/id%d.dat.sz", root);
               data_bits = readfrombinary_char(binfile_sz, &data_bytes);		
            }

            MPI_Bcast(&data_bytes, 1, MPI_INT, root, MPI_COMM_WORLD);
            compress_ratio += data_bytes*8.0/(size*sizeof(double)*8);
         
            if(id != root)
            {
                  data_bits = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes);
            }
            MPI_Bcast(data_bits, data_bytes, MPI_UNSIGNED_CHAR, root, MPI_COMM_WORLD);

            char binfile_zs[64];
            sprintf(binfile_zs, "dataset/id%d.dat.zs", root);
            writetobinary_char(binfile_zs, data_bits, data_bytes); //.dat.zs
            char sz_decomp_cmd[64];
            sprintf(sz_decomp_cmd, "%sdataset/id%d%s%d", sz_decomp_cmd_prefix, root, sz_decomp_cmd_suffix, size);
            //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
            int iret_decomp = system(sz_decomp_cmd); //.dat.zs --> .dat.zs.out
            char binfile_out[64];
            sprintf(binfile_out, "dataset/id%d.dat.zs.out", root);
            char txtfile[64];
            sprintf(txtfile, "dataset/id%d.dat.zs.out.txt", root); 
            double* decompressed_data = readfrombinary_writetotxt_double(binfile_out, txtfile, size);	

            double gs = 0;
            for(int i=0; i<size; i++)
            {
                  if(id == root)
                  {
                     gs += fabs(decompressed_data[i] - save[i]);
                  }
                  else
                  {
                     save[i] = decompressed_data[i];
                  }
            } 
            gosa = gs/size;    

            //todo
            free(data_bits);
         }		
         else if(CT == 1)
         {
            int array_double_len;
            struct vector msg; 

            if(id == root)
            {
                  //mycommpress
                  double* array_double = NULL;
                  char* array_char = NULL;
                  int* array_char_displacement = NULL;
                  array_double_len = myCompress_double(save, &array_double, &array_char, &array_char_displacement, size);
                  msg.p_data = array_double;
                  msg.c_data = array_char;
                  msg.disp = array_char_displacement;
            }

            MPI_Bcast(&array_double_len, 1, MPI_INT, root, MPI_COMM_WORLD);
            int num_p = array_double_len, num_c = size - array_double_len;
            compress_ratio += (float)(num_c*sizeof(char)+num_p*sizeof(double))/((num_c+num_p)*sizeof(double));
         
            if(id != root)
            {
                  msg.p_data = (double*) malloc(sizeof(double)*num_p);
                  if(num_c > 0)
                  {
                     msg.c_data = (char*) malloc(sizeof(char)*num_c);
                     msg.disp = (int*) malloc(sizeof(int)*num_c);					
                  }
            }
            MPI_Bcast(msg.p_data, num_p, MPI_DOUBLE, root, MPI_COMM_WORLD);
            if(num_c == 0)
            {
                  msg.c_data = (char*) malloc(sizeof(char)*1);
                  msg.disp = (int*) malloc(sizeof(int)*1);
                  msg.c_data[0] = 'z';
                  msg.disp[0] = -1;
            }	
            else
            {
                  MPI_Bcast(msg.c_data, num_c, MPI_CHAR, root, MPI_COMM_WORLD);
                  MPI_Bcast(msg.disp, num_c, MPI_INT, root, MPI_COMM_WORLD);	
            }

            double* decompressed_data = myDecompress_double(msg.p_data, msg.c_data, msg.disp, size);
            double gs = 0;
            for(int i=0; i<size; i++)
            {
                  if(id == root)
                  {
                     gs += fabs(decompressed_data[i]-save[i]);
                  }
                  else
                  {
                     save[i] = decompressed_data[i];
                  }
            }
            gosa += gs/size;

            //todo
            free(msg.p_data);
            free(msg.c_data);
            free(msg.disp);
         }
         else if(CT == 0)
         {
            MPI_Bcast(save, mx_size - diag_ref, MPI_DOUBLE, j % p, MPI_COMM_WORLD);
         }         
      }
   }

   double end = MPI_Wtime();

   if (id == root_p) {
      // printf("[LU]\n");
      // print_mx(A, mx_size * mx_size, mx_size);
      printf("\n[L]\n");
      L_print(A, mx_size);
      printf("\n[U]\n");
      U_print(A, mx_size);
      ln();
      printf("mpi: %f s\n", end - start);
      ln();

		printf("--------------------------------------------------\n");
		printf("FINAL RESULTS:\n");	

		//printf("id = %d, elapsed = %f = %f - %f\n", id, end_time-start_time, end_time, start_time);
      int loop = mx_size*(mx_size-1)/2;
		printf("gosa = %f \n", gosa/loop);
		printf("compression ratio: sz %f, nolossy_performance %f, nolossy_area %f \n", 1/(sz_comp_ratio/loop), 1/(nolossy_performance/loop), 1/(nolossy_area/loop));
		printf("compress ratio = %f \n", 1/(compress_ratio/loop));   
      printf("resent = %d (percentage = %f)\n", resent, 1.0*resent/((p-1)*loop));  
      //printf("p=%d, tmp_size=%d, mx_size=%d, diag_ref=%d\n", p, tmp_size, mx_size, diag_ref); 

      char fn[] = "lu.csv";
      int fexist = access(fn, 0);
      FILE* fp = fopen(fn, "a"); 
      if(fexist == -1)
      {
         fprintf(fp, "nprocs, matrix size, CT, absErrorBound, BER, compression ratio, time, gosa, resent, resent ratio\n"); 
      }    
      fprintf(fp, "%d, %d, %d, %e, %e, %f, %f, %f, %d, %f\n", p, mx_size, CT, absErrorBound, BER, 1/(compress_ratio/loop), end - start, gosa/loop, resent, 1.0*resent/((p-1)*loop));    
      fclose(fp);    
   }
   free(A);

   MPI_Finalize();
   return EXIT_SUCCESS;
}

/*
 * gen_mx - generate contiguous matrix
 *
 * @dim dim x dim matrix
 * @return matrix
 */
double *gen_mx (size_t dim)
{
   int i, j, tot = dim * dim;
   double *M = malloc(sizeof(double) * tot);
   for (i = 0; i < tot; i++) {
      M[i] = rand() % 101 - 50;
   }

   return M;
}

/*
 * mx_print - dumb matrix print function
 *
 * @M matrix/row
 * @dim matrix/row dimension
 * @sep where put separator
 */
void print_mx (double *M, size_t dim, size_t sep)
{
   int i, j;
   for (i = 0; i < dim; i++) {
      printf("% *.*f\t", 4, 2, M[i]);
      if ((i + 1) % sep == 0) {
         ln();
      }
   }
}

/*
 * forw_elim - forward Gauss elimination
 *
 * @origin row pointer by reference
 * @master_row row in which lays diagonal
 */
void forw_elim(double **origin, double *master_row, size_t dim)
{
   if (**origin == 0)
      return;

   double k = **origin / master_row[0];

   int i;
   for (i = 1; i < dim; i++) {
      (*origin)[i] = (*origin)[i] - k * master_row[i];
   }
   **origin = k;
}

/*
 * U_print - dumb U matrix print function
 */
void U_print (double *M, int dim)
{
   int i, j;
   double z = 0;
   for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
         if (j >= i) {
            printf("% *.*f\t", 4, 2, M[i * dim + j]);
         } else {
            printf("% *.*f\t", 4, 2, z);
         }
      }
      ln();
   }
}

/*
 * L_print - dumb L matrix print function
 */
void L_print (double *M, int dim)
{
   int i, j;
   double z = 0, u = 1;
   for (i = 0; i < dim; i++) {
      for (j = 0; j < dim; j++) {
         if (j > i) {
            printf("% *.*f\t", 4, 2, z);
         } else if (i == j) {
            printf("% *.*f\t", 4, 2, u);
         } else {
            printf("% *.*f\t", 4, 2, M[i * dim + j]);
         }
      }
      ln();
   }
}
