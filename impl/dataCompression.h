/*
 *
 */
#define absErrBound         0.1 //default 0.0001=2^{-12} (-13?), 0.000001=2^{-20}, 0.00001=2^{-16}, 0.001=2^{-10}, 0.01=2^{-7}
#define absErrBound_binary  20 //SZ, equal to above
#define relBoundRatio       0.01
#define pw_relBoundRatio    0.01    
#define CT                  1 //compress type for pingpong & himeno & k-means, 0 no compress, 1 mycompress, 2 no-lossy-performance, 3 no-lossy-area, 4 sz
#define byte_or_bit         1 //1 byte, 2 bit
#define data_num            8192 //pingpong
#define filename            "testdouble_8_8_8_128.txt" //k-means, "input.txt", "testdouble_8_8_128.txt", "testdouble_8_8_8_128.txt"

float calcCompressionRatio_himeno_ij_ik_jk(float[MIMAX][MJMAX][MKMAX], int, int);
// MPI_Datatype myCompress_himeno(void*, int, int, int, int, int, int);
float calcCompressionRatio_himeno_sz(float[MIMAX][MJMAX][MKMAX], int, int);
float calcCompressionRatio_himeno_nolossy_performance(float[MIMAX][MJMAX][MKMAX], int, int);
float calcCompressionRatio_himeno_nolossy_area(float[MIMAX][MJMAX][MKMAX], int, int);
void getFloatBin(float, char[]);
float* readFileFloat(char[]);
int myCompress(float[], float**, char**, int**, int);
float* myDecompress(float[], char[], int[], int);
int myCompress_double(double[], double**, char**, int**, int);
double* myDecompress_double(double[], char[], int[], int);
float* transform_3d_array_to_1d_array(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int, int);