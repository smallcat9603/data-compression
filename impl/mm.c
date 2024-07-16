#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <malloc.h>
#include <mpi.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include "param.h"
#include "dataCompression.h"

//todo
struct vector
{
  double* p_data; //precise data
  char* c_data; //compressed data
  int* disp; //displacement of compressed data
};

typedef struct {
        unsigned int rows;
        unsigned int cols;
        double **mat_data;
} matrix_struct;

matrix_struct *get_matrix_struct(char matrix[]) {
    matrix_struct *m = malloc(sizeof(matrix_struct));
    m->rows = 0;
    m->cols = 0;
    FILE* myfile = fopen(matrix, "r");
    
    if(myfile == NULL) {
        printf("Error: The file you entered could not be found.\n");
        exit(EXIT_FAILURE);
    }
    // get the rows and columns
    int ch = 0;
    do {
        ch = fgetc(myfile);
        
        // count the columns at the first line (looking for "\t")
        if(m->rows == 0 && ch == '\t')
            m->cols++;
        
        // count the rows with "\n"
        if(ch == '\n')
            m->rows++;
            
    } while (ch != EOF);
    
    // write rows and cols to struct
    m->cols++;
    
    // allocate memory for matrix data
    m->mat_data = calloc(m->rows, sizeof(double*)); 
    int i;
    for(i=0; i < m->rows; ++i)
        m->mat_data[i]=calloc(m->cols, sizeof(double));
        
    
    rewind(myfile);
    int x,y;
    
    // fill matrix with data
    for(x = 0; x < m->rows; x++) {
        for(y = 0; y < m->cols; y++) {
            if (!fscanf(myfile, "%lf", &m->mat_data[x][y])) 
            break;
        }
    }
    
    fclose(myfile);

    return m;
}

void print_matrix(matrix_struct *matrix_to_print){
    int i,j;
    for(i = 0; i < matrix_to_print->rows; i++) {
        for(j = 0; j < matrix_to_print->cols; j++) {
            printf("%lf\t",matrix_to_print->mat_data[i][j]); //Use lf format specifier, \n is for new line
        }
        printf("\n");
    }
}

void free_matrix(matrix_struct *matrix_to_free) {
    for(int i = 0; i < matrix_to_free->rows; i++) {
        free(matrix_to_free->mat_data[i]);
    }
    free(matrix_to_free->mat_data);
    free(matrix_to_free);
}

double *mat_2D_to_1D(matrix_struct *m) {
    double *matrix = malloc( (m->rows * m->cols) * sizeof(double) );
    for (int i = 0; i < m->rows; i++) {
        memcpy( matrix + (i * m->cols), m->mat_data[i], m->cols * sizeof(double) );
    }
    return matrix;
}

int main(int argc, char *argv[]) {
	//modify CT
	int CT = 0;
	if(argc > 3) CT = atoi(argv[3]);

    /** Matrix Properties
     * [0] = Rows of Matrix A
     * [1] = Cols of Matrix A
     * [2] = Rows of Matrix B
     * [3] = Cols of Matrix B
     **/
    int matrix_properties[4];
     
    double *m_a = NULL;
    double *m_b = NULL;
    double *final_matrix = NULL;
    
    int num_worker, rank;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_worker);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    /** the master initializes the data **/
    if (rank == 0) {
        
        if(argc < 3){
            printf("ERROR: Please specify only 2 files.\n");
            exit(EXIT_FAILURE);
        }
            
        matrix_struct *m_1 = get_matrix_struct(argv[1]);
        matrix_struct *m_2 = get_matrix_struct(argv[2]);

        if(m_1->cols != m_2->rows){
            printf("ERROR: The number of columns of matrix A must match the number of rows of matrix B.\n");
            exit(EXIT_FAILURE);
        }
        
        if (m_1->rows % num_worker != 0) {
            printf("ERROR: Matrix can not be calculated with this number of tasks.\n");
            exit(EXIT_FAILURE);
        }
        
        // fill the property-array for workers
        matrix_properties[0] = m_1->rows;
        matrix_properties[1] = m_1->cols;
        matrix_properties[2] = m_2->rows;
        matrix_properties[3] = m_2->cols;
        
        /* generate 1D matrices for workers 
         * m_a is the 1D Matrix of m_1 
         * m_a is the 1D Matrix of m_1 
        */
        m_a = mat_2D_to_1D(m_1);
        m_b = mat_2D_to_1D(m_2);

        free_matrix(m_1);
        free_matrix(m_2);
    }

    // send the matrix properties to the workers
    MPI_Bcast(&matrix_properties, 4, MPI_INT, 0, MPI_COMM_WORLD);

    // calculate the 1D-sizes of the matrices
    int size_a   = matrix_properties[0] * matrix_properties[1];
    int size_b   = matrix_properties[2] * matrix_properties[3];
    int size_res = matrix_properties[0] * matrix_properties[3];
    
    // allocate memory for 1D-matrices
    if(rank == 0) {
        final_matrix = malloc( size_res * sizeof(double) );
    } else {
        m_a = malloc( size_a * sizeof(double) );
        m_b = malloc( size_b * sizeof(double) );
    }

	double gosa = 0;
	double start_time, end_time;
	float compress_ratio = 0;
	float sz_comp_ratio = 0;
	float nolossy_performance = 0;
	float nolossy_area = 0;
    uint32_t crc_a = 0;
    uint32_t crc_check_a = 0;
    unsigned char crc_ok_a = 'y';
    unsigned char* crc_ok_recv_a = NULL;
    uint32_t crc_b = 0;
    uint32_t crc_check_b = 0;
    unsigned char crc_ok_b = 'y';
    unsigned char* crc_ok_recv_b = NULL;    
    int resent = 0;
    srand((unsigned)time(NULL));      
	start_time = MPI_Wtime();
    

    if(CT == 10)
    {
        MPI_Bcast_bitwise_crc_hamming(m_a, size_a, 0, rank, num_worker, &compress_ratio, &gosa, &resent);
        MPI_Bcast_bitwise_crc_hamming(m_b, size_b, 0, rank, num_worker, &compress_ratio, &gosa, &resent);  
    }
    else if(CT == 9)
    {
        MPI_Bcast_bitwise_mask_crc(m_a, size_a, 0, rank, num_worker, &compress_ratio, &gosa, &resent);
        MPI_Bcast_bitwise_mask_crc(m_b, size_b, 0, rank, num_worker, &compress_ratio, &gosa, &resent);
    }
    else if(CT == 8)
    {
        MPI_Bcast_bitwise_crc(m_a, size_a, 0, rank, num_worker, &compress_ratio, &gosa, &resent);
        MPI_Bcast_bitwise_crc(m_b, size_b, 0, rank, num_worker, &compress_ratio, &gosa, &resent);  

        // int data_bytes_a = 0, data_bytes_b = 0;
        // double a_min = 0, b_min = 0;

        // unsigned char* data_bits_a = NULL;
        // unsigned char* data_bits_b = NULL;
        
        // //a
        // if(rank == 0)
        // {
        //     // sz_comp_ratio += calcCompressionRatio_sz_double(m_a, size_a);
        //     // nolossy_performance += calcCompressionRatio_nolossy_performance_double(m_a, size_a);
        //     // nolossy_area += calcCompressionRatio_nolossy_area_double(m_a, size_a);

        //     //mycommpress
        //     double* a_small = NULL;
        //     a_min = toSmallDataset_double(m_a, &a_small, size_a);

        //     int data_pos_a = 8; //position of filled bit in last byte --> 87654321

        //     myCompress_bitwise_double(a_small, size_a, &data_bits_a, &data_bytes_a, &data_pos_a);	
        //     crc_a = do_crc32(data_bits_a, data_bytes_a);		
        // }

        // MPI_Bcast(&data_bytes_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Bcast(&a_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // compress_ratio += data_bytes_a*8.0/(size_a*sizeof(double)*8);
    
        // if(rank != 0)
        // {
        //     data_bits_a = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_a);
        // }
        // MPI_Bcast(data_bits_a, data_bytes_a, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        // MPI_Bcast(&crc_a, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        // if(rank != 0)
        // {
        //     crc_check_a = do_crc32(data_bits_a, data_bytes_a);

        //     if(BER > 0)
        //     {
        //         double ber = BER;
        //         uint64_t to = 1/ber;
        //         uint64_t r = get_random_int(0, to);
        //         if(r < data_bytes_a * 8)
        //         {
        //             crc_check_a = 0;
        //         }
        //     }
            
        //     if (crc_a == crc_check_a)
        //     {
        //         // printf("CRC passed\n");
        //         crc_ok_a = 'y';
        //     }  
        //     else
        //     {
        //         // printf("CRC NOT passed\n");
        //         crc_ok_a = 'n';
        //     }            
        // }
        // else
        // {
        //     crc_ok_recv_a = (unsigned char *)malloc(num_worker*1*sizeof(unsigned char));
        // }

        // MPI_Gather(&crc_ok_a, 1, MPI_UNSIGNED_CHAR, crc_ok_recv_a, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        
        // if(rank == 0)
        // {
        //     for(int i = 1; i < num_worker; i++)
        //     {
        //         if(crc_ok_recv_a[i] == 'n')
        //         {
        //             MPI_Send(data_bits_a, data_bytes_a, MPI_UNSIGNED_CHAR, i, i, MPI_COMM_WORLD);
        //             resent++;
        //         }
        //     }
        // }
        // else if(crc_ok_a == 'n')
        // {
        //     MPI_Recv(data_bits_a, data_bytes_a, MPI_UNSIGNED_CHAR, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // }

        // //b
        // if(rank == 0)
        // {
        //     // sz_comp_ratio += calcCompressionRatio_sz_double(m_b, size_b);
        //     // nolossy_performance += calcCompressionRatio_nolossy_performance_double(m_b, size_b);
        //     // nolossy_area += calcCompressionRatio_nolossy_area_double(m_b, size_b);

        //     //mycommpress
        //     double* b_small = NULL;
        //     b_min = toSmallDataset_double(m_b, &b_small, size_b);

        //     int data_pos_b = 8; //position of filled bit in last byte --> 87654321

        //     myCompress_bitwise_double(b_small, size_b, &data_bits_b, &data_bytes_b, &data_pos_b);	
        //     crc_b = do_crc32(data_bits_b, data_bytes_b);		
        // }

        // MPI_Bcast(&data_bytes_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        // MPI_Bcast(&b_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        // compress_ratio += data_bytes_b*8.0/(size_b*sizeof(double)*8);
    
        // if(rank != 0)
        // {
        //     data_bits_b = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_b);
        // }
        // MPI_Bcast(data_bits_b, data_bytes_b, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        // MPI_Bcast(&crc_b, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

        // if(rank != 0)
        // {
        //     crc_check_b = do_crc32(data_bits_b, data_bytes_b);

        //     if(BER > 0)
        //     {
        //         double ber = BER;
        //         uint64_t to = 1/ber;
        //         uint64_t r = get_random_int(0, to);
        //         if(r < data_bytes_b * 8)
        //         {
        //             crc_check_b = 0;
        //         }
        //     }
            
        //     if (crc_b == crc_check_b)
        //     {
        //         // printf("CRC passed\n");
        //         crc_ok_b = 'y';
        //     }  
        //     else
        //     {
        //         // printf("CRC NOT passed\n");
        //         crc_ok_b = 'n';
        //     }            
        // }
        // else
        // {
        //     crc_ok_recv_b = (unsigned char *)malloc(num_worker*1*sizeof(unsigned char));
        // }

        // MPI_Gather(&crc_ok_b, 1, MPI_UNSIGNED_CHAR, crc_ok_recv_b, 1, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);
        
        // if(rank == 0)
        // {
        //     for(int i = 1; i < num_worker; i++)
        //     {
        //         if(crc_ok_recv_b[i] == 'n')
        //         {
        //             MPI_Send(data_bits_b, data_bytes_b, MPI_UNSIGNED_CHAR, i, i, MPI_COMM_WORLD);
        //             resent++;
        //         }
        //     }
        // }
        // else if(crc_ok_b == 'n')
        // {
        //     MPI_Recv(data_bits_b, data_bytes_b, MPI_UNSIGNED_CHAR, 0, rank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        // }

        // double* decompressed_data_a = myDecompress_bitwise_double(data_bits_a, data_bytes_a, size_a);
        // double* decompressed_data_b = myDecompress_bitwise_double(data_bits_b, data_bytes_b, size_b);
        // for(int i=0; i<size_a; i++)
        // {
        //     if(rank == 0)
        //     {
        //         gosa += fabs(decompressed_data_a[i] + a_min - m_a[i]);
        //     }
        //     else
        //     {
        //         m_a[i] = decompressed_data_a[i] + a_min;
        //     }
        // }
        // for(int i=0; i<size_b; i++)
        // {
        //     if(rank == 0)
        //     {
        //         gosa += fabs(decompressed_data_b[i] + b_min - m_b[i]);
        //     }
        //     else
        //     {
        //         m_b[i] = decompressed_data_b[i] + b_min;
        //     }
        // }

        // //todo
        // free(data_bits_a);
        // free(data_bits_b);
    }    
    else if(CT == 7)
    {
        int data_bytes_a = 0, data_bytes_b = 0;
        double a_min = 0, b_min = 0;

        unsigned char* data_bits_a = NULL;
        unsigned char* data_bits_b = NULL;

        int type_a = 0, type_b = 0;
        double medium_a = 0, medium_b = 0;
        
        //a
        if(rank == 0)
        {
            //mycommpress
            double* a_small = NULL;
            a_min = toSmallDataset_double(m_a, &a_small, size_a);

            int data_pos_a = 8; //position of filled bit in last byte --> 87654321

            medium_a = med_dataset_double(a_small, size_a, &type_a);
            char double_arr[64+1];
            doubletostr(&medium_a, double_arr);
            char mask[1+11+8];
            strncpy(mask, double_arr, 1+11+8);			

            myCompress_bitwise_double_mask(a_small, size_a, &data_bits_a, &data_bytes_a, &data_pos_a, type_a, mask);			
        }

        MPI_Bcast(&data_bytes_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&a_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_a*8.0/(size_a*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_a = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_a);
        }
        MPI_Bcast(data_bits_a, data_bytes_a, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        MPI_Bcast(&medium_a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&type_a, 1, MPI_INT, 0, MPI_COMM_WORLD);

        //b
        if(rank == 0)
        {
            //mycommpress
            double* b_small = NULL;
            b_min = toSmallDataset_double(m_b, &b_small, size_b);

            int data_pos_b = 8; //position of filled bit in last byte --> 87654321

            medium_b = med_dataset_double(b_small, size_b, &type_b);
            char double_arr[64+1];
            doubletostr(&medium_b, double_arr);
            char mask[1+11+8];
            strncpy(mask, double_arr, 1+11+8);					

            myCompress_bitwise_double_mask(b_small, size_b, &data_bits_b, &data_bytes_b, &data_pos_b, type_b, mask);				
        }

        MPI_Bcast(&data_bytes_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&b_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_b*8.0/(size_b*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_b = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_b);
        }
        MPI_Bcast(data_bits_b, data_bytes_b, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        MPI_Bcast(&medium_b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&type_b, 1, MPI_INT, 0, MPI_COMM_WORLD);

        char double_arr_a_recv[64+1];
        char double_arr_b_recv[64+1];
        doubletostr(&medium_a, double_arr_a_recv);
        doubletostr(&medium_b, double_arr_b_recv);
        char mask_a_recv[1+11+8];
        char mask_b_recv[1+11+8];
        strncpy(mask_a_recv, double_arr_a_recv, 1+11+8);
        strncpy(mask_b_recv, double_arr_b_recv, 1+11+8);    			

        double* decompressed_data_a = myDecompress_bitwise_double_mask(data_bits_a, data_bytes_a, size_a, type_a, mask_a_recv);
        double* decompressed_data_b = myDecompress_bitwise_double_mask(data_bits_b, data_bytes_b, size_b, type_b, mask_b_recv);
        double gs = 0;
        for(int i=0; i<size_a; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_a[i] + a_min - m_a[i]);
            }
            else
            {
                m_a[i] = decompressed_data_a[i] + a_min;
            }
        }
        gosa += gs/size_a;
        gs = 0;
        for(int i=0; i<size_b; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_b[i] + b_min - m_b[i]);
            }
            else
            {
                m_b[i] = decompressed_data_b[i] + b_min;
            }
        }
        gosa += gs/size_b;        

        //todo
        free(data_bits_a);
        free(data_bits_b);
    }
    else if(CT == 6)
    {
        int data_bytes_a = 0, data_bytes_b = 0;
        double a_min = 0, b_min = 0;

        unsigned char* data_bits_a = NULL;
        unsigned char* data_bits_b = NULL;
        
        //a
        if(rank == 0)
        {
            //mycommpress
            double* a_small = NULL;
            a_min = toSmallDataset_double(m_a, &a_small, size_a);

            int data_pos_a = 8; //position of filled bit in last byte --> 87654321

            myCompress_bitwise_double_np(a_small, size_a, &data_bits_a, &data_bytes_a, &data_pos_a);			
        }

        MPI_Bcast(&data_bytes_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&a_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_a*8.0/(size_a*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_a = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_a);
        }
        MPI_Bcast(data_bits_a, data_bytes_a, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        //b
        if(rank == 0)
        {
            //mycommpress
            double* b_small = NULL;
            b_min = toSmallDataset_double(m_b, &b_small, size_b);

            int data_pos_b = 8; //position of filled bit in last byte --> 87654321

            myCompress_bitwise_double_np(b_small, size_b, &data_bits_b, &data_bytes_b, &data_pos_b);			
        }

        MPI_Bcast(&data_bytes_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&b_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_b*8.0/(size_b*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_b = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_b);
        }
        MPI_Bcast(data_bits_b, data_bytes_b, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        double* decompressed_data_a = myDecompress_bitwise_double_np(data_bits_a, data_bytes_a, size_a);
        double* decompressed_data_b = myDecompress_bitwise_double_np(data_bits_b, data_bytes_b, size_b);
        double gs = 0;
        for(int i=0; i<size_a; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_a[i] + a_min - m_a[i]);
            }
            else
            {
                m_a[i] = decompressed_data_a[i] + a_min;
            }
        }
        gosa += gs/size_a; 
        gs = 0;
        for(int i=0; i<size_b; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_b[i] + b_min - m_b[i]);
            }
            else
            {
                m_b[i] = decompressed_data_b[i] + b_min;
            }
        }
        gosa += gs/size_b;          

        //todo
        free(data_bits_a);
        free(data_bits_b);
    }
    else if(CT == 5)
    {
        int data_bytes_a = 0, data_bytes_b = 0;
        double a_min = 0, b_min = 0;

        unsigned char* data_bits_a = NULL;
        unsigned char* data_bits_b = NULL;
        
        //a
        if(rank == 0)
        {
            // sz_comp_ratio += calcCompressionRatio_sz_double(m_a, size_a);
            // nolossy_performance += calcCompressionRatio_nolossy_performance_double(m_a, size_a);
            // nolossy_area += calcCompressionRatio_nolossy_area_double(m_a, size_a);

            //mycommpress
            double* a_small = NULL;
            a_min = toSmallDataset_double(m_a, &a_small, size_a);

            int data_pos_a = 8; //position of filled bit in last byte --> 87654321

            myCompress_bitwise_double(a_small, size_a, &data_bits_a, &data_bytes_a, &data_pos_a);			
        }

        MPI_Bcast(&data_bytes_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&a_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_a*8.0/(size_a*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_a = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_a);
        }
        MPI_Bcast(data_bits_a, data_bytes_a, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        //b
        if(rank == 0)
        {
            // sz_comp_ratio += calcCompressionRatio_sz_double(m_b, size_b);
            // nolossy_performance += calcCompressionRatio_nolossy_performance_double(m_b, size_b);
            // nolossy_area += calcCompressionRatio_nolossy_area_double(m_b, size_b);

            //mycommpress
            double* b_small = NULL;
            b_min = toSmallDataset_double(m_b, &b_small, size_b);

            int data_pos_b = 8; //position of filled bit in last byte --> 87654321

            myCompress_bitwise_double(b_small, size_b, &data_bits_b, &data_bytes_b, &data_pos_b);			
        }

        MPI_Bcast(&data_bytes_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&b_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_b*8.0/(size_b*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_b = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_b);
        }
        MPI_Bcast(data_bits_b, data_bytes_b, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        double* decompressed_data_a = myDecompress_bitwise_double(data_bits_a, data_bytes_a, size_a);
        double* decompressed_data_b = myDecompress_bitwise_double(data_bits_b, data_bytes_b, size_b);
        double gs = 0;
        for(int i=0; i<size_a; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_a[i] + a_min - m_a[i]);
            }
            else
            {
                m_a[i] = decompressed_data_a[i] + a_min;
            }
        }
        gosa += gs/size_a; 
        gs = 0;
        for(int i=0; i<size_b; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_b[i] + b_min - m_b[i]);
            }
            else
            {
                m_b[i] = decompressed_data_b[i] + b_min;
            }
        }
        gosa += gs/size_b; 

        //todo
        free(data_bits_a);
        free(data_bits_b);
    }
    else if(CT == 4)
    {
        int data_bytes_a = 0, data_bytes_b = 0;

        unsigned char* data_bits_a = NULL;
        unsigned char* data_bits_b = NULL;
        
        //a
        if(rank == 0)
        {
            char binfile[64];
            sprintf(binfile, "dataset/a.dat");
            writetobinary_double(binfile, m_a, size_a); //.txt --> .dat
            char sz_comp_cmd[64];
            sprintf(sz_comp_cmd, "%s%g%sdataset/a%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, sz_comp_cmd_suffix2, size_a);
            //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
            int iret_comp = system(sz_comp_cmd); //.dat --> .dat.sz
            char binfile_sz[64];
            sprintf(binfile_sz, "dataset/a.dat.sz");
            data_bits_a = readfrombinary_char(binfile_sz, &data_bytes_a);		
        }

        MPI_Bcast(&data_bytes_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_a*8.0/(size_a*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_a = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_a);
        }
        MPI_Bcast(data_bits_a, data_bytes_a, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        //b
        if(rank == 0)
        {
            char binfile[64];
            sprintf(binfile, "dataset/b.dat");
            writetobinary_double(binfile, m_b, size_b); //.txt --> .dat
            char sz_comp_cmd[64];
            sprintf(sz_comp_cmd, "%s%g%sdataset/b%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, sz_comp_cmd_suffix2, size_b);
            //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
            int iret_comp = system(sz_comp_cmd); //.dat --> .dat.sz
            char binfile_sz[64];
            sprintf(binfile_sz, "dataset/b.dat.sz");
            data_bits_b = readfrombinary_char(binfile_sz, &data_bytes_b);				
        }

        MPI_Bcast(&data_bytes_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        compress_ratio += data_bytes_b*8.0/(size_b*sizeof(double)*8);
    
        if(rank != 0)
        {
            data_bits_b = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_b);
        }
        MPI_Bcast(data_bits_b, data_bytes_b, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

        char binfile_zs_a[64];
        sprintf(binfile_zs_a, "dataset/a.dat.zs");
        writetobinary_char(binfile_zs_a, data_bits_a, data_bytes_a); //.dat.zs
        char sz_decomp_cmd_a[64];
        sprintf(sz_decomp_cmd_a, "%sdataset/a%s%d", sz_decomp_cmd_prefix, sz_decomp_cmd_suffix, size_a);
        //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
        int iret_decomp_a = system(sz_decomp_cmd_a); //.dat.zs --> .dat.zs.out
        char binfile_out_a[64];
        sprintf(binfile_out_a, "dataset/a.dat.zs.out");
        char txtfile_a[64];
        sprintf(txtfile_a, "dataset/a.dat.zs.out.txt"); 
        double* decompressed_data_a = readfrombinary_writetotxt_double(binfile_out_a, txtfile_a, size_a);			

        char binfile_zs_b[64];
        sprintf(binfile_zs_b, "dataset/b.dat.zs");
        writetobinary_char(binfile_zs_b, data_bits_b, data_bytes_b); //.dat.zs
        char sz_decomp_cmd_b[64];
        sprintf(sz_decomp_cmd_b, "%sdataset/b%s%d", sz_decomp_cmd_prefix, sz_decomp_cmd_suffix, size_b);
        //int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
        int iret_decomp_b = system(sz_decomp_cmd_b); //.dat.zs --> .dat.zs.out
        char binfile_out_b[64];
        sprintf(binfile_out_b, "dataset/b.dat.zs.out");
        char txtfile_b[64];
        sprintf(txtfile_b, "dataset/b.dat.zs.out.txt");  
        double* decompressed_data_b = readfrombinary_writetotxt_double(binfile_out_b, txtfile_b, size_b);	

        double gs = 0;
        for(int i=0; i<size_a; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_a[i] - m_a[i]);
            }
            else
            {
                m_a[i] = decompressed_data_a[i];
            }
        }
        gosa += gs/size_a; 
        gs = 0;
        for(int i=0; i<size_b; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_b[i] - m_b[i]);
            }
            else
            {
                m_b[i] = decompressed_data_b[i];
            }
        }
        gosa += gs/size_b;          

        //todo
        free(data_bits_a);
        free(data_bits_b);
    }		
    else if(CT == 1)
    {
        int array_double_len_a, array_double_len_b;
        struct vector msg_a, msg_b; 

        //a
        if(rank == 0)
        {
            //mycommpress
            double* array_double_a = NULL;
            char* array_char_a = NULL;
            int* array_char_displacement_a = NULL;
            array_double_len_a = myCompress_double(m_a, &array_double_a, &array_char_a, &array_char_displacement_a, size_a);
            msg_a.p_data = array_double_a;
            msg_a.c_data = array_char_a;
            msg_a.disp = array_char_displacement_a;
        }

        MPI_Bcast(&array_double_len_a, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int num_p_a = array_double_len_a, num_c_a = size_a - array_double_len_a;
        compress_ratio += (float)(num_c_a*sizeof(char)+num_p_a*sizeof(double))/((num_c_a+num_p_a)*sizeof(double));
    
        if(rank != 0)
        {
            msg_a.p_data = (double*) malloc(sizeof(double)*num_p_a);
            if(num_c_a > 0)
            {
                msg_a.c_data = (char*) malloc(sizeof(char)*num_c_a);
                msg_a.disp = (int*) malloc(sizeof(int)*num_c_a);					
            }
        }
        MPI_Bcast(msg_a.p_data, num_p_a, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(num_c_a == 0)
        {
            msg_a.c_data = (char*) malloc(sizeof(char)*1);
            msg_a.disp = (int*) malloc(sizeof(int)*1);
            msg_a.c_data[0] = 'z';
            msg_a.disp[0] = -1;
        }	
        else
        {
            MPI_Bcast(msg_a.c_data, num_c_a, MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Bcast(msg_a.disp, num_c_a, MPI_INT, 0, MPI_COMM_WORLD);	
        }
        
        //b
        if(rank == 0)
        {
            double* array_double_b = NULL;
            char* array_char_b = NULL;
            int* array_char_displacement_b = NULL;
            array_double_len_b = myCompress_double(m_b, &array_double_b, &array_char_b, &array_char_displacement_b, size_b);
            msg_b.p_data = array_double_b;
            msg_b.c_data = array_char_b;
            msg_b.disp = array_char_displacement_b;
        }

        MPI_Bcast(&array_double_len_b, 1, MPI_INT, 0, MPI_COMM_WORLD);
        int num_p_b = array_double_len_b, num_c_b = size_b - array_double_len_b;
        compress_ratio += (float)(num_c_b*sizeof(char)+num_p_b*sizeof(double))/((num_c_b+num_p_b)*sizeof(double));

        if(rank != 0)
        {
            msg_b.p_data = (double*) malloc(sizeof(double)*num_p_b);
            if(num_c_b > 0)
            {
                msg_b.c_data = (char*) malloc(sizeof(char)*num_c_b);
                msg_b.disp = (int*) malloc(sizeof(int)*num_c_b);
            }
        }					
        MPI_Bcast(msg_b.p_data, num_p_b, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if(num_c_b == 0)
        {
            msg_b.c_data = (char*) malloc(sizeof(char)*1);
            msg_b.disp = (int*) malloc(sizeof(int)*1);
            msg_b.c_data[0] = 'z';
            msg_b.disp[0] = -1;
        }	
        else
        {			
            MPI_Bcast(msg_b.c_data, num_c_b, MPI_CHAR, 0, MPI_COMM_WORLD);
            MPI_Bcast(msg_b.disp, num_c_b, MPI_INT, 0, MPI_COMM_WORLD);
        }

        double gs = 0; 

        double* decompressed_data_a = myDecompress_double(msg_a.p_data, msg_a.c_data, msg_a.disp, size_a);
        for(int i=0; i<size_a; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_a[i]-m_a[i]);
            }
            else
            {
                m_a[i] = decompressed_data_a[i];
            }
        }
        gosa += gs/size_a; 
        gs = 0;
        double* decompressed_data_b = myDecompress_double(msg_b.p_data, msg_b.c_data, msg_b.disp, size_b);
        for(int i=0; i<size_b; i++)
        {
            if(rank == 0)
            {
                gs += fabs(decompressed_data_b[i]-m_b[i]);
            }
            else
            {
                m_b[i] = decompressed_data_b[i];
            }				
        }
        gosa += gs/size_b; 

        //todo
        free(msg_a.p_data);
        free(msg_a.c_data);
        free(msg_a.disp);
        free(msg_b.p_data);
        free(msg_b.c_data);
        free(msg_b.disp);
    }
    else if(CT == 0)
    {
        // send 1D matrices to workers
        MPI_Bcast(m_a, size_a , MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(m_b, size_b , MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    end_time = MPI_Wtime();
    
    // calculate the start- and endrow for worker  
    int startrow = rank * ( matrix_properties[0] / num_worker);
    int endrow = ((rank + 1) * ( matrix_properties[0] / num_worker)) -1;
    
    /* calculate sub matrices */
    int number_of_rows = size_res / num_worker;
    double *result_matrix = calloc(number_of_rows, sizeof(double));

    int position = 0;

    for (int i = startrow; i <= endrow; i++) {
        for (int j = 0; j < matrix_properties[3]; j++) {
            for (int k = 0; k < matrix_properties[2]; k++) {
                result_matrix[position] +=
                    m_a[ (i * matrix_properties[1] + k) ] *
                    m_b[ (k * matrix_properties[3] + j) ];
            }
            position++;
        }
    }
    
    free(m_a);
    free(m_b);
    
    /* collect the results */
    MPI_Gather(result_matrix, number_of_rows, MPI_DOUBLE,
           final_matrix, number_of_rows,  MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /** The master presents the results on the console */
    if (rank == 0){
        // int size = matrix_properties[0] * matrix_properties[3];
        // int i = 0;
        // while (i < size) {
        //     printf("%lf\t", final_matrix[i]);
        //     i++;
        
        //     if (i % matrix_properties[3] == 0)
        //         printf("\n");
        // }

		printf("--------------------------------------------------\n");
		printf("FINAL RESULTS:\n");	

		printf("rank = %d, elapsed = %f = %f - %f\n", rank, end_time-start_time, end_time, start_time);
		printf("gosa = %f \n", gosa/2);
		printf("compression ratio: sz %f, nolossy_performance %f, nolossy_area %f \n", 1/(sz_comp_ratio/2), 1/(nolossy_performance/2), 1/(nolossy_area/2));
		printf("compress ratio = %f \n", 1/(compress_ratio/2));    
        printf("resent = %d (percentage = %f)\n", resent, resent/(2.0*(num_worker-1)));    
        // printf("hamming_correct = %d (percentage = %f)\n", hamming_correct, 1.0*hamming_correct/(resent+hamming_correct));

        char fn[] = "mm.csv";
        int fexist = access(fn, 0);
        FILE* fp = fopen(fn, "a"); 
        if(fexist == -1)
        {
            fprintf(fp, "num_worker, size_res, CT, absErrorBound, BER, compression ratio, time, gosa, resent, resent ratio\n"); 
        }    
        fprintf(fp, "%d, %d, %d, %e, %e, %f, %f, %f, %d, %f\n", num_worker, size_res, CT, absErrorBound, BER, 1/(compress_ratio/2), end_time - start_time, gosa/2, resent, resent/(2.0*(num_worker-1)));    
        fclose(fp);          
    }
    
    free(result_matrix);
    free(final_matrix);
    
    MPI_Finalize();
    exit(EXIT_SUCCESS);
}
