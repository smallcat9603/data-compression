/*
 *
 */
#define absErrorBound	0.000001 //default 0.0001=2^{-12} (-13?), 0.000001=2^{-20}, 0.00001=2^{-16}, 0.001=2^{-10}, 0.01=2^{-7}
// #define absErrorBound_binary  20 //bitwise, SZ, equal to above
/*compress type 
0 no compress, 
1 mycompress, 
2 no-lossy-performance, 
3 no-lossy-area, 
4 sz, 
5 bitwise mycompress, 
6 bitwise no prediction, 
7 bitmask-based bitwise, 
8 bitwise w/ crc,
9 bitmask-based bitwise w/crc  
10 bitwise w/ crc hamming 
11 bitwise only prediction 
*/

int MPI_Bcast_bitwise_double(void *buf, int count, MPI_Datatype datatype, int root, MPI_Comm comm);

int MPI_Send_bitwise_double(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Recv_bitwise_double(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Send_bitwise_double_np(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Recv_bitwise_double_np(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);
int MPI_Send_bitwise_double_op(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm);
int MPI_Recv_bitwise_double_op(void *buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status);

double* myDecompress_bitwise_double_mask(unsigned char*, int, int, int, char[1+11+8]);
double decompress_bitwise_double_mask(char*, int, double, double, double, int, char[1+11+8]);
void compress_bitwise_double_mask(double, unsigned char**, int*, int*, int, char[1+11+8]);
void myCompress_bitwise_double_mask(double[], int, unsigned char**, int*, int*, int, char[1+11+8]);
float* myDecompress_bitwise_mask(unsigned char*, int, int, int, char[1+8+8]);
float decompress_bitwise_float_mask(char*, int, float, float, float, int, char[1+8+8]);
void compress_bitwise_float_mask(float, unsigned char**, int*, int*, int, char[1+8+8]);
void myCompress_bitwise_mask(float[], int, unsigned char**, int*, int*, int, char[1+8+8]);

double* myDecompress_bitwise_double_np(unsigned char*, int, int);
double decompress_bitwise_double_np(char*, int);
float* myDecompress_bitwise_np(unsigned char*, int, int);
float decompress_bitwise_float_np(char*, int);
void myCompress_bitwise_double_np(double[], int, unsigned char**, int*, int*);
void myCompress_bitwise_np(float[], int, unsigned char**, int*, int*);

void myCompress_bitwise_double_op(double[], int, unsigned char**, int*, int*);
double* myDecompress_bitwise_double_op(unsigned char*, int, int);
void myCompress_bitwise_op(float[], int, unsigned char**, int*, int*);
float* myDecompress_bitwise_op(unsigned char*, int, int);

double* myDecompress_bitwise_double(unsigned char*, int, int);
double decompress_bitwise_double(char*, int, double, double, double);
float* myDecompress_bitwise(unsigned char*, int, int);
float decompress_bitwise_float(char*, int, float, float, float);

void myCompress_bitwise_double(double[], int, unsigned char**, int*, int*);
void myCompress_bitwise(float[], int, unsigned char**, int*, int*);

void compress_bitwise_double(double, unsigned char**, int*, int*);
void compress_bitwise_float(float, unsigned char**, int*, int*);

double toSmallDataset_double(double[], double**, int);
float toSmallDataset_float(float[], float**, int);

double med_dataset_double(double*, int, int*);
float med_dataset_float(float*, int, int*);

void getFloatBin(float, char[]);
void getDoubleBin(double,char[]);
float* readFileFloat(char[]);
int myCompress(float[], float**, char**, int**, int);
float* myDecompress(float[], char[], int[], int);
int myCompress_double(double[], double**, char**, int**, int);
double* myDecompress_double(double[], char[], int[], int);
void floattostr(float*, char*);
void doubletostr(double*, char*);
float strtofloat(char*);
double strtodbl(char*);
void writetobinary_float(const char*, float*, int);
void writetobinary_double(const char*, double*, int);
void writetobinary_char(const char *, unsigned char*, int);
float* readfrombinary_float(const char*, int);
double* readfrombinary_double(const char*, int);
unsigned char* readfrombinary_char(const char *, int*);
float* readfrombinary_writetotxt_float(const char*, const char*, int);
double* readfrombinary_writetotxt_double(const char*, const char*, int);

void add_bit_to_bytes(unsigned char**, int*, int*, int);
void bit_set(unsigned char*, unsigned char, int);

int to_absErrorBound_binary(double absErrBound);

void cast_bits_to_char(unsigned char* bits, char* data, int bytes); //1 bit --> char '0' or '1' (8 bits)
