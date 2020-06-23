/*  This is an implementation of the k-means clustering algorithm (aka Lloyd's algorithm) using MPI (message passing interface). */
//Original version: https://github.com/dzdao/k-means-clustering-mpi

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<unistd.h>
#include<math.h>
#include<errno.h>
#include<mpi.h>
#include<string.h>
#include "param.h"
#include "dataCompression.h"

#define MAX_ITERATIONS 1000

int numOfClusters = 0;
int numOfElements = 0;
int num_of_processes = 0;
// numOfElements%num_of_processes should be 0 for MPI_Scatter

//todo
struct vector
{
  double* p_data; //precise data
  char* c_data; //compressed data
  int* disp; //displacement of compressed data
};

/* This function goes through that data points and assigns them to a cluster */
void assign2Cluster(double k_x[], double k_y[], double recv_x[], double recv_y[], int assign[])
{
	double min_dist = 10000000;
	double x=0, y=0, temp_dist=0;
	int k_min_index = 0;

	for(int i = 0; i < (numOfElements/num_of_processes)/* + 1*/; i++)
	{
		//fix bug
		x = fabs(recv_x[i] - k_x[0]);
		y = fabs(recv_y[i] - k_y[0]);
		min_dist = sqrt((x*x) + (y*y));
		k_min_index = 0;

		for(int j = 1; j < numOfClusters; j++)
		{
			x = fabs(recv_x[i] - k_x[j]);
			y = fabs(recv_y[i] - k_y[j]);
			temp_dist = sqrt((x*x) + (y*y));

			// new minimum distance found
			if(temp_dist < min_dist)
			{
				min_dist = temp_dist;
				k_min_index = j;
			}
		}

		// update the cluster assignment of this data points
		assign[i] = k_min_index;
	}

}

/* Recalcuate k-means of each cluster because each data point may have
   been reassigned to a new cluster for each iteration of the algorithm */
void calcKmeans(double k_means_x[], double k_means_y[], double data_x_points[], double data_y_points[], int k_assignment[])
{
	double total_x = 0;
	double total_y = 0;
	int numOfpoints = 0;

	for(int i = 0; i < numOfClusters; i++)
	{
		total_x = 0;
		total_y = 0;
		numOfpoints = 0;

		for(int j = 0; j < numOfElements; j++)
		{
			if(k_assignment[j] == i)
			{
				total_x += data_x_points[j];
				total_y += data_y_points[j];
				numOfpoints++;
			}
		}

		if(numOfpoints != 0)
		{
			k_means_x[i] = total_x / numOfpoints;
			k_means_y[i] = total_y / numOfpoints;
		}
	}

}

int main(int argc, char *argv[])
{
	// initialize the MPI environment
	MPI_Init(NULL, NULL);

	// get number of processes
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	// get rank
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	//printf("my rank: %d \n", world_rank);

	// send buffers
	double *k_means_x = NULL;		// k means corresponding x values
	double *k_means_y = NULL;		// k means corresponding y values
	int *k_assignment = NULL;		// each data point is assigned to a cluster
	double *data_x_points = NULL;
	double *data_y_points = NULL;

	// receive buffer
	double *recv_x = NULL;
	double *recv_y = NULL;
	int *recv_assign = NULL;

	if(world_rank == 0)
	{
		if(argc != 2)
		{
			printf("Please include an argument after the program name to list how many processes.\n");
			printf("e.g. To indicate 4 processes, run: mpirun -n 4 ./kmeans 4\n");
			exit(-1);
		}

		num_of_processes = atoi(argv[1]);

		// char buffer[2];
		// printf("How many clusters would you like to analyze for? ");
		// scanf("%s", buffer);
		// printf("\n");

		// numOfClusters = atoi(buffer);
		numOfClusters = clusters;
		printf("Ok %d clusters it is.\n", numOfClusters);

		// broadcast the number of clusters to all nodes
		MPI_Bcast(&numOfClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// allocate memory for arrays
		k_means_x = (double *)malloc(sizeof(double) * numOfClusters);
		k_means_y = (double *)malloc(sizeof(double) * numOfClusters);

		if(k_means_x == NULL || k_means_y == NULL)
		{
			perror("malloc");
			exit(-1);
		}

		printf("Reading input data from file...\n\n");

		char* filename_suffix = filename suffix;
		FILE* fp = fopen(filename_suffix, "r");

		if(!fp)
		{
			perror("fopen");
			exit(-1);
		}

		// count number of lines to find out how many elements
		int c = 0;
		numOfElements = 0;
		while(!feof(fp))
		{
			c = fgetc(fp);
			if(c == '\n')
			{
				numOfElements++;
			}
		}

		//x, y
		numOfElements = numOfElements/2;

		printf("There are a total number of %d elements in the file.\n", numOfElements);

		// broadcast the number of elements to all nodes
		MPI_Bcast(&numOfElements, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// allocate memory for an array of data points
		data_x_points = (double *)malloc(sizeof(double) * numOfElements);
		data_y_points = (double *)malloc(sizeof(double) * numOfElements);
		k_assignment = (int *)malloc(sizeof(int) * numOfElements);

		if(data_x_points == NULL || data_y_points == NULL || k_assignment == NULL)
		{
			perror("malloc");
			exit(-1);
		}

		// reset file pointer to origin of file
		fseek(fp, 0, SEEK_SET);

		// now read in points and fill the arrays
		int i = 0;

		double point_x=0, point_y=0;

		while(fscanf(fp, "%lf %lf", &point_x, &point_y) != EOF)
		{
			data_x_points[i] = point_x;
			data_y_points[i] = point_y;

			// assign the initial k means to zero
			k_assignment[i] = 0;
			i++;
		}

		// close file pointer
		fclose(fp);

		// randomly select initial k-means
		time_t t;
		srand((unsigned) time(&t));
		int random;
		for(int i = 0; i < numOfClusters; i++) {
			//random = rand() % numOfElements;
			random = numOfElements/numOfClusters * i;
			k_means_x[i] = data_x_points[random];
			k_means_y[i] = data_y_points[random];
		}

		printf("Running k-means algorithm for %d iterations...\n\n", MAX_ITERATIONS);
		// for(int i = 0; i < numOfClusters; i++)
		// {
		// 	printf("Initial K-means: (%f, %f)\n", k_means_x[i], k_means_y[i]);
		// }

		// allocate memory for receive buffers
		recv_x = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes)/* + 1*/));
		recv_y = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes)/* + 1*/));
		recv_assign = (int *)malloc(sizeof(int) * ((numOfElements/num_of_processes)/* + 1*/));

		if(recv_x == NULL || recv_y == NULL || recv_assign == NULL)
		{
			perror("malloc");
			exit(-1);
		}
	}
	else
	{	// I am a worker node

		num_of_processes = atoi(argv[1]);

		// receive broadcast of number of clusters
		MPI_Bcast(&numOfClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// receive broadcast of number of elements
		MPI_Bcast(&numOfElements, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// allocate memory for arrays
		k_means_x = (double *)malloc(sizeof(double) * numOfClusters);
		k_means_y = (double *)malloc(sizeof(double) * numOfClusters);

		if(k_means_x == NULL || k_means_y == NULL)
		{
			perror("malloc");
			exit(-1);
		}

		// allocate memory for receive buffers
		recv_x = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes)/* + 1*/));
		recv_y = (double *)malloc(sizeof(double) * ((numOfElements/num_of_processes)/* + 1*/));
		recv_assign = (int *)malloc(sizeof(int) * ((numOfElements/num_of_processes)/* + 1*/));

		if(recv_x == NULL || recv_y == NULL || recv_assign == NULL)
		{
			perror("malloc");
			exit(-1);
		}
	}

	/* Distribute the work among all nodes. The data points itself will stay constant and
	   not change for the duration of the algorithm. */
	MPI_Scatter(data_x_points, (numOfElements/num_of_processes)/* + 1*/, MPI_DOUBLE,
		recv_x, (numOfElements/num_of_processes)/* + 1*/, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	MPI_Scatter(data_y_points, (numOfElements/num_of_processes)/* + 1*/, MPI_DOUBLE,
		recv_y, (numOfElements/num_of_processes)/* + 1*/, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	int count = 0;
	double gosa = 0;
	double start_time, end_time;
	float compress_ratio = 0;
	float sz_comp_ratio = 0;
	float nolossy_performance = 0;
	float nolossy_area = 0;
	start_time = MPI_Wtime();
	while(count < MAX_ITERATIONS)
	{
		if(CT == 7)
		{
			int data_bytes_x = 0, data_bytes_y = 0;
			double k_means_x_min = 0, k_means_y_min = 0;

			unsigned char* data_bits_x = NULL;
			unsigned char* data_bits_y = NULL;

			int type_x = 0, type_y = 0;
			double medium_x = 0, medium_y = 0;
			
			//x
			if(world_rank == 0)
			{
				//mycommpress
				double* k_means_x_small = NULL;
				k_means_x_min = toSmallDataset_double(k_means_x, &k_means_x_small, numOfClusters);

				int data_pos_x = 8; //position of filled bit in last byte --> 87654321

				medium_x = med_dataset_double(k_means_x_small, numOfClusters, &type_x);
				char double_arr[64+1];
				doubletostr(&medium_x, double_arr);
				char mask[1+11+8];
				strncpy(mask, double_arr, 1+11+8);			

				myCompress_bitwise_double_mask(k_means_x_small, numOfClusters, &data_bits_x, &data_bytes_x, &data_pos_x, type_x, mask);			
			}

			MPI_Bcast(&data_bytes_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&k_means_x_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_x*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_x = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_x);
			}
			MPI_Bcast(data_bits_x, data_bytes_x, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			MPI_Bcast(&medium_x, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&type_x, 1, MPI_INT, 0, MPI_COMM_WORLD);

			//y
			if(world_rank == 0)
			{
				//mycommpress
				double* k_means_y_small = NULL;
				k_means_y_min = toSmallDataset_double(k_means_y, &k_means_y_small, numOfClusters);

				int data_pos_y = 8; //position of filled bit in last byte --> 87654321

				medium_y = med_dataset_double(k_means_y_small, numOfClusters, &type_y);
				char double_arr[64+1];
				doubletostr(&medium_y, double_arr);
				char mask[1+11+8];
				strncpy(mask, double_arr, 1+11+8);					

				myCompress_bitwise_double_mask(k_means_y_small, numOfClusters, &data_bits_y, &data_bytes_y, &data_pos_y, type_y, mask);				
			}

			MPI_Bcast(&data_bytes_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&k_means_y_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_y*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_y = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_y);
			}
			MPI_Bcast(data_bits_y, data_bytes_y, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			MPI_Bcast(&medium_y, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&type_y, 1, MPI_INT, 0, MPI_COMM_WORLD);

			char double_arr_x_recv[64+1];
			char double_arr_y_recv[64+1];
			doubletostr(&medium_x, double_arr_x_recv);
			doubletostr(&medium_y, double_arr_y_recv);
			char mask_x_recv[1+11+8];
			char mask_y_recv[1+11+8];
			strncpy(mask_x_recv, double_arr_x_recv, 1+11+8);
			strncpy(mask_y_recv, double_arr_y_recv, 1+11+8);    			

			double* decompressed_data_x = myDecompress_bitwise_double_mask(data_bits_x, data_bytes_x, numOfClusters, type_x, mask_x_recv);
			double* decompressed_data_y = myDecompress_bitwise_double_mask(data_bits_y, data_bytes_y, numOfClusters, type_y, mask_y_recv);
			for(int i=0; i<numOfClusters; i++)
			{
				if(world_rank == 0)
				{
					gosa += fabs(decompressed_data_x[i] + k_means_x_min - k_means_x[i]);
					gosa += fabs(decompressed_data_y[i] + k_means_y_min - k_means_y[i]);
				}
				else
				{
					k_means_x[i] = decompressed_data_x[i] + k_means_x_min;
					k_means_y[i] = decompressed_data_y[i] + k_means_y_min;
				}
			}
			
			if(count == MAX_ITERATIONS - 1)
			{
				gosa = gosa/(2*MAX_ITERATIONS*numOfClusters);
			}

			//todo
			free(data_bits_x);
			free(data_bits_y);
		}

		if(CT == 6)
		{
			int data_bytes_x = 0, data_bytes_y = 0;
			double k_means_x_min = 0, k_means_y_min = 0;

			unsigned char* data_bits_x = NULL;
			unsigned char* data_bits_y = NULL;
			
			//x
			if(world_rank == 0)
			{
				// sz_comp_ratio += calcCompressionRatio_sz_double(k_means_x, numOfClusters);
				// nolossy_performance += calcCompressionRatio_nolossy_performance_double(k_means_x, numOfClusters);
				// nolossy_area += calcCompressionRatio_nolossy_area_double(k_means_x, numOfClusters);

				//mycommpress
				double* k_means_x_small = NULL;
				k_means_x_min = toSmallDataset_double(k_means_x, &k_means_x_small, numOfClusters);

				int data_pos_x = 8; //position of filled bit in last byte --> 87654321

				myCompress_bitwise_double_np(k_means_x_small, numOfClusters, &data_bits_x, &data_bytes_x, &data_pos_x);			
			}

			MPI_Bcast(&data_bytes_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&k_means_x_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_x*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_x = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_x);
			}
			MPI_Bcast(data_bits_x, data_bytes_x, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			//y
			if(world_rank == 0)
			{
				// sz_comp_ratio += calcCompressionRatio_sz_double(k_means_y, numOfClusters);
				// nolossy_performance += calcCompressionRatio_nolossy_performance_double(k_means_y, numOfClusters);
				// nolossy_area += calcCompressionRatio_nolossy_area_double(k_means_y, numOfClusters);

				//mycommpress
				double* k_means_y_small = NULL;
				k_means_y_min = toSmallDataset_double(k_means_y, &k_means_y_small, numOfClusters);

				int data_pos_y = 8; //position of filled bit in last byte --> 87654321

				myCompress_bitwise_double_np(k_means_y_small, numOfClusters, &data_bits_y, &data_bytes_y, &data_pos_y);			
			}

			MPI_Bcast(&data_bytes_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&k_means_y_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_y*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_y = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_y);
			}
			MPI_Bcast(data_bits_y, data_bytes_y, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			double* decompressed_data_x = myDecompress_bitwise_double_np(data_bits_x, data_bytes_x, numOfClusters);
			double* decompressed_data_y = myDecompress_bitwise_double_np(data_bits_y, data_bytes_y, numOfClusters);
			for(int i=0; i<numOfClusters; i++)
			{
				if(world_rank == 0)
				{
					gosa += fabs(decompressed_data_x[i] + k_means_x_min - k_means_x[i]);
					gosa += fabs(decompressed_data_y[i] + k_means_y_min - k_means_y[i]);
				}
				else
				{
					k_means_x[i] = decompressed_data_x[i] + k_means_x_min;
					k_means_y[i] = decompressed_data_y[i] + k_means_y_min;
				}
			}
			
			if(count == MAX_ITERATIONS - 1)
			{
				gosa = gosa/(2*MAX_ITERATIONS*numOfClusters);
			}

			//todo
			free(data_bits_x);
			free(data_bits_y);
		}

		if(CT == 5)
		{
			int data_bytes_x = 0, data_bytes_y = 0;
			double k_means_x_min = 0, k_means_y_min = 0;

			unsigned char* data_bits_x = NULL;
			unsigned char* data_bits_y = NULL;
			
			//x
			if(world_rank == 0)
			{
				// sz_comp_ratio += calcCompressionRatio_sz_double(k_means_x, numOfClusters);
				// nolossy_performance += calcCompressionRatio_nolossy_performance_double(k_means_x, numOfClusters);
				// nolossy_area += calcCompressionRatio_nolossy_area_double(k_means_x, numOfClusters);

				//mycommpress
				double* k_means_x_small = NULL;
				k_means_x_min = toSmallDataset_double(k_means_x, &k_means_x_small, numOfClusters);

				int data_pos_x = 8; //position of filled bit in last byte --> 87654321

				myCompress_bitwise_double(k_means_x_small, numOfClusters, &data_bits_x, &data_bytes_x, &data_pos_x);			
			}

			MPI_Bcast(&data_bytes_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&k_means_x_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_x*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_x = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_x);
			}
			MPI_Bcast(data_bits_x, data_bytes_x, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			//y
			if(world_rank == 0)
			{
				// sz_comp_ratio += calcCompressionRatio_sz_double(k_means_y, numOfClusters);
				// nolossy_performance += calcCompressionRatio_nolossy_performance_double(k_means_y, numOfClusters);
				// nolossy_area += calcCompressionRatio_nolossy_area_double(k_means_y, numOfClusters);

				//mycommpress
				double* k_means_y_small = NULL;
				k_means_y_min = toSmallDataset_double(k_means_y, &k_means_y_small, numOfClusters);

				int data_pos_y = 8; //position of filled bit in last byte --> 87654321

				myCompress_bitwise_double(k_means_y_small, numOfClusters, &data_bits_y, &data_bytes_y, &data_pos_y);			
			}

			MPI_Bcast(&data_bytes_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
			MPI_Bcast(&k_means_y_min, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_y*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_y = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_y);
			}
			MPI_Bcast(data_bits_y, data_bytes_y, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			double* decompressed_data_x = myDecompress_bitwise_double(data_bits_x, data_bytes_x, numOfClusters);
			double* decompressed_data_y = myDecompress_bitwise_double(data_bits_y, data_bytes_y, numOfClusters);
			for(int i=0; i<numOfClusters; i++)
			{
				if(world_rank == 0)
				{
					gosa += fabs(decompressed_data_x[i] + k_means_x_min - k_means_x[i]);
					gosa += fabs(decompressed_data_y[i] + k_means_y_min - k_means_y[i]);
				}
				else
				{
					k_means_x[i] = decompressed_data_x[i] + k_means_x_min;
					k_means_y[i] = decompressed_data_y[i] + k_means_y_min;
				}
			}
			
			if(count == MAX_ITERATIONS - 1)
			{
				gosa = gosa/(2*MAX_ITERATIONS*numOfClusters);
			}

			//todo
			free(data_bits_x);
			free(data_bits_y);
		}

		if(CT == 4)
		{
			int data_bytes_x = 0, data_bytes_y = 0;

			unsigned char* data_bits_x = NULL;
			unsigned char* data_bits_y = NULL;
			
			//x
			if(world_rank == 0)
			{
				char binfile[64];
    			sprintf(binfile, "dataset/x%d.dat", count);
				writetobinary_double(binfile, k_means_x, numOfClusters); //.txt --> .dat
				char sz_comp_cmd[64];
				sprintf(sz_comp_cmd, "%s%g%sdataset/x%d%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, count, sz_comp_cmd_suffix2, numOfClusters);
				//int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
				int iret_comp = system(sz_comp_cmd); //.dat --> .dat.sz
				char binfile_sz[64];
				sprintf(binfile_sz, "dataset/x%d.dat.sz", count);
				data_bits_x = readfrombinary_char(binfile_sz, &data_bytes_x);		
			}

			MPI_Bcast(&data_bytes_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_x*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_x = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_x);
			}
			MPI_Bcast(data_bits_x, data_bytes_x, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			//y
			if(world_rank == 0)
			{
				char binfile[64];
				sprintf(binfile, "dataset/y%d.dat", count);
				writetobinary_double(binfile, k_means_y, numOfClusters); //.txt --> .dat
				char sz_comp_cmd[64];
				sprintf(sz_comp_cmd, "%s%g%sdataset/y%d%s%d", sz_comp_cmd_prefix, absErrorBound, sz_comp_cmd_suffix1, count, sz_comp_cmd_suffix2, numOfClusters);
				//int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
				int iret_comp = system(sz_comp_cmd); //.dat --> .dat.sz
				char binfile_sz[64];
				sprintf(binfile_sz, "dataset/y%d.dat.sz", count);
				data_bits_y = readfrombinary_char(binfile_sz, &data_bytes_y);				
			}

			MPI_Bcast(&data_bytes_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
			compress_ratio += data_bytes_y*8.0/(numOfClusters*sizeof(double)*8);
		
			if(world_rank != 0)
			{
				data_bits_y = (unsigned char*) malloc(sizeof(unsigned char)*data_bytes_y);
			}
			MPI_Bcast(data_bits_y, data_bytes_y, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

			char binfile_zs_x[64];
			sprintf(binfile_zs_x, "dataset/x%d.dat.zs", count);
			writetobinary_char(binfile_zs_x, data_bits_x, data_bytes_x); //.dat.zs
			char sz_decomp_cmd_x[64];
			sprintf(sz_decomp_cmd_x, "%sdataset/x%d%s%d", sz_decomp_cmd_prefix, count, sz_decomp_cmd_suffix, numOfClusters);
			//int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
			int iret_decomp_x = system(sz_decomp_cmd_x); //.dat.zs --> .dat.zs.out
			char binfile_out_x[64];
			sprintf(binfile_out_x, "dataset/x%d.dat.zs.out", count);
			char txtfile_x[64];
			sprintf(txtfile_x, "dataset/x%d.dat.zs.out.txt", count); 
			double* decompressed_data_x = readfrombinary_writetotxt_double(binfile_out_x, txtfile_x, numOfClusters);			

			char binfile_zs_y[64];
			sprintf(binfile_zs_y, "dataset/y%d.dat.zs", count);
			writetobinary_char(binfile_zs_y, data_bits_y, data_bytes_y); //.dat.zs
			char sz_decomp_cmd_y[64];
			sprintf(sz_decomp_cmd_y, "%sdataset/y%d%s%d", sz_decomp_cmd_prefix, count, sz_decomp_cmd_suffix, numOfClusters);
			//int iret = system("./sz -z -f -c sz.config -M ABS -A 0.001 -i ./testdata/x86/testfloat_8_8_128.dat -1 8192");
			int iret_decomp_y = system(sz_decomp_cmd_y); //.dat.zs --> .dat.zs.out
			char binfile_out_y[64];
			sprintf(binfile_out_y, "dataset/y%d.dat.zs.out", count);
			char txtfile_y[64];
			sprintf(txtfile_y, "dataset/y%d.dat.zs.out.txt", count);  
			double* decompressed_data_y = readfrombinary_writetotxt_double(binfile_out_y, txtfile_y, numOfClusters);	

			for(int i=0; i<numOfClusters; i++)
			{
				if(world_rank == 0)
				{
					gosa += fabs(decompressed_data_x[i] - k_means_x[i]);
					gosa += fabs(decompressed_data_y[i] - k_means_y[i]);
				}
				else
				{
					k_means_x[i] = decompressed_data_x[i];
					k_means_y[i] = decompressed_data_y[i];
				}
			}
			
			if(count == MAX_ITERATIONS - 1)
			{
				gosa = gosa/(2*MAX_ITERATIONS*numOfClusters);
			}

			//todo
			free(data_bits_x);
			free(data_bits_y);
		}		

		if(CT == 1)
		{
			int array_double_len_x, array_double_len_y;
			struct vector msg_x, msg_y; 

			//x
			if(world_rank == 0)
			{
				//mycommpress
				double* array_double_x = NULL;
				char* array_char_x = NULL;
				int* array_char_displacement_x = NULL;
				array_double_len_x = myCompress_double(k_means_x, &array_double_x, &array_char_x, &array_char_displacement_x, numOfClusters);
				msg_x.p_data = array_double_x;
				msg_x.c_data = array_char_x;
				msg_x.disp = array_char_displacement_x;
			}

			MPI_Bcast(&array_double_len_x, 1, MPI_INT, 0, MPI_COMM_WORLD);
			int num_p_x = array_double_len_x, num_c_x = numOfClusters - array_double_len_x;
			compress_ratio += (float)(num_c_x*sizeof(char)+num_p_x*sizeof(double))/((num_c_x+num_p_x)*sizeof(double));
		
			if(world_rank != 0)
			{
				msg_x.p_data = (double*) malloc(sizeof(double)*num_p_x);
				if(num_c_x > 0)
				{
					msg_x.c_data = (char*) malloc(sizeof(char)*num_c_x);
					msg_x.disp = (int*) malloc(sizeof(int)*num_c_x);					
				}
			}
			MPI_Bcast(msg_x.p_data, num_p_x, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if(num_c_x == 0)
			{
				msg_x.c_data = (char*) malloc(sizeof(char)*1);
				msg_x.disp = (int*) malloc(sizeof(int)*1);
				msg_x.c_data[0] = 'z';
				msg_x.disp[0] = -1;
			}	
			else
			{
				MPI_Bcast(msg_x.c_data, num_c_x, MPI_CHAR, 0, MPI_COMM_WORLD);
				MPI_Bcast(msg_x.disp, num_c_x, MPI_INT, 0, MPI_COMM_WORLD);	
			}
			
			//y
			if(world_rank == 0)
			{
				double* array_double_y = NULL;
				char* array_char_y = NULL;
				int* array_char_displacement_y = NULL;
				array_double_len_y = myCompress_double(k_means_y, &array_double_y, &array_char_y, &array_char_displacement_y, numOfClusters);
				msg_y.p_data = array_double_y;
				msg_y.c_data = array_char_y;
				msg_y.disp = array_char_displacement_y;
			}

			MPI_Bcast(&array_double_len_y, 1, MPI_INT, 0, MPI_COMM_WORLD);
			int num_p_y = array_double_len_y, num_c_y = numOfClusters - array_double_len_y;
			compress_ratio += (float)(num_c_y*sizeof(char)+num_p_y*sizeof(double))/((num_c_y+num_p_y)*sizeof(double));

			if(world_rank != 0)
			{
				msg_y.p_data = (double*) malloc(sizeof(double)*num_p_y);
				if(num_c_y > 0)
				{
					msg_y.c_data = (char*) malloc(sizeof(char)*num_c_y);
					msg_y.disp = (int*) malloc(sizeof(int)*num_c_y);
				}
			}					
			MPI_Bcast(msg_y.p_data, num_p_y, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			if(num_c_y == 0)
			{
				msg_y.c_data = (char*) malloc(sizeof(char)*1);
				msg_y.disp = (int*) malloc(sizeof(int)*1);
				msg_y.c_data[0] = 'z';
				msg_y.disp[0] = -1;
			}	
			else
			{			
				MPI_Bcast(msg_y.c_data, num_c_y, MPI_CHAR, 0, MPI_COMM_WORLD);
				MPI_Bcast(msg_y.disp, num_c_y, MPI_INT, 0, MPI_COMM_WORLD);
			}

			double* decompressed_data_x = myDecompress_double(msg_x.p_data, msg_x.c_data, msg_x.disp, numOfClusters);
			for(int i=0; i<numOfClusters; i++)
			{
				if(world_rank == 0)
				{
					gosa += fabs(decompressed_data_x[i]-k_means_x[i]);
				}
				else
				{
					k_means_x[i] = decompressed_data_x[i];
				}
			}

			double* decompressed_data_y = myDecompress_double(msg_y.p_data, msg_y.c_data, msg_y.disp, numOfClusters);
			for(int i=0; i<numOfClusters; i++)
			{
				if(world_rank == 0)
				{
					gosa += fabs(decompressed_data_y[i]-k_means_y[i]);
				}
				else
				{
					k_means_y[i] = decompressed_data_y[i];
				}				
			}
			
			if(count == MAX_ITERATIONS - 1)
			{
				gosa = gosa/(2*MAX_ITERATIONS*numOfClusters);
			}

			//todo
			free(msg_x.p_data);
			free(msg_x.c_data);
			free(msg_x.disp);
			free(msg_y.p_data);
			free(msg_y.c_data);
			free(msg_y.disp);
		}		

		if(CT == 0)
		{
			// broadcast k-means arrays
			MPI_Bcast(k_means_x, numOfClusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Bcast(k_means_y, numOfClusters, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		// scatter k-cluster assignments array
		MPI_Scatter(k_assignment, (numOfElements/num_of_processes)/* + 1*/, MPI_INT,
			recv_assign, (numOfElements/num_of_processes)/* + 1*/, MPI_INT, 0, MPI_COMM_WORLD);

		// assign the data points to a cluster
		assign2Cluster(k_means_x, k_means_y, recv_x, recv_y, recv_assign);

		// gather back k-cluster assignments
		MPI_Gather(recv_assign, (numOfElements/num_of_processes)/* + 1*/, MPI_INT,
			k_assignment, (numOfElements/num_of_processes)/* + 1*/, MPI_INT, 0, MPI_COMM_WORLD);

		// let the root process recalculate k means
		if(world_rank == 0)
		{
			calcKmeans(k_means_x, k_means_y, data_x_points, data_y_points, k_assignment);
			//printf("Finished iteration %d\n",count);
		}

		count++;
	}
	end_time = MPI_Wtime();

	if(world_rank == 0)
	{
		printf("--------------------------------------------------\n");
		printf("FINAL RESULTS:\n");
		// for(int i = 0; i < numOfClusters; i++)
		// {
		// 	printf("Cluster #%d: (%f, %f)\n", i, k_means_x[i], k_means_y[i]);
		// }
		// printf("--------------------------------------------------\n");

		// for(int i = 0; i < numOfElements; i++)
		// {
		// 	printf("%d, ", k_assignment[i]);
		// }
		// printf("\n--------------------------------------------------\n");

		//char* output_filename = filename output_suffix;
		char output_filename[32];
		sprintf(output_filename, "%s%s%d_%g_%d%s", filename, output_suffix, CT, absErrorBound, clusters, suffix);
		FILE* fp = fopen(output_filename, "w");
		for(int i = 0; i < numOfElements; i++)
		{
			fprintf(fp, "%d\n", k_assignment[i]);
		}
		for(int i = 0; i < numOfClusters; i++)
		{
			//printf("Cluster #%d: (%f, %f)\n", i, k_means_x[i], k_means_y[i]);
			fprintf(fp, "%lf\n", k_means_x[i]);
			fprintf(fp, "%lf\n", k_means_y[i]);

		}
		printf("--------------------------------------------------\n");
		fclose(fp);		

		printf("rank = %d, elapsed = %f = %f - %f\n", world_rank, end_time-start_time, end_time, start_time);
		printf("gosa = %f \n", gosa);
		printf("compression ratio: sz %f, nolossy_performance %f, nolossy_area %f \n", 1/(sz_comp_ratio/(2*MAX_ITERATIONS)), 1/(nolossy_performance/(2*MAX_ITERATIONS)), 1/(nolossy_area/(2*MAX_ITERATIONS)));
		printf("compress ratio = %f \n", 1/(compress_ratio/(2*MAX_ITERATIONS)));
	}

	// deallocate memory and clean up
	free(k_means_x);
	free(k_means_y);
	free(data_x_points);
	free(data_y_points);
	free(k_assignment);
	free(recv_x);
	free(recv_y);
	free(recv_assign);

	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

}
