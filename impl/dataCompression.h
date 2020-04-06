/*
 *
 */
#define absErrBound         0.000001 //default 0.0001=2^{-12} (-13?), 0.000001=2^{-20}, 0.00001=2^{-16}, 0.001=2^{-10}, 0.01=2^{-7}
#define absErrBound_binary  20 //bitwise, SZ, equal to above
#define relBoundRatio       0.01
#define pw_relBoundRatio    0.01    
#define CT                  1 //compress type for pingpong & himeno & k-means, 0 no compress, 1 mycompress, 2 no-lossy-performance, 3 no-lossy-area, 4 sz
#define byte_or_bit         1 //1 byte, 2 bit
#define data_num            8192 //pingpong
#define filename            "dataset/testfloat_8_8_128" //pingpong, k-means, "input", "testfloat_8_8_128", "testdouble_8_8_128", "testdouble_8_8_8_128", obs_info, num_plasma
#define suffix              ".txt" //k-means
#define output_suffix       "_output_" //k-means
#define clusters            100 //k-means

void myCompress_bitwise(float[], int, unsigned char**, int&, int&);

void compress_bitwise_double(double, unsigned char**, int&, int&);
void compress_bitwise_float(float, unsigned char**, int&, int&);

double toSmallDataset(double[], double**, int);
float toSmallDataset(float[], float**, int);

float calCompressRatio_bitwise_float(float[], int);
float calCompressRatio_bitwise_double(double[], int);
float calCompressRatio_bitwise_double2(float[], int);

float calcCompressionRatio_himeno_ij_ik_jk(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
// MPI_Datatype myCompress_himeno(void*, int, int, int, int, int, int);
float calcCompressionRatio_himeno_sz(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
float calcCompressionRatio_himeno_nolossy_performance(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
float calcCompressionRatio_himeno_nolossy_area(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
void getFloatBin(float, char[]);
float* readFileFloat(char[]);
int myCompress(float[], float**, char**, int**, int);
float* myDecompress(float[], char[], int[], int);
int myCompress_double(double[], double**, char**, int**, int);
double* myDecompress_double(double[], char[], int[], int);
float* transform_3d_array_to_1d_array(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
void floattostr(float*, char*);
void doubletostr(double*, char*);
float strtofloat(char*);
double strtodbl(char*);
void writetobinary_float(const char*, float*, int);
void writetobinary_double(const char*, double*, int);
float* readfrombinary_float(const char*, int);
double* readfrombinary_double(const char*, int);
void readfrombinary_writetotxt_float(const char*, const char*, int);
void readfrombinary_writetotxt_double(const char*, const char*, int);

void add_bit_to_bytes(unsigned char**, int&, int&, int);
void bit_set(unsigned char*, unsigned char, int);