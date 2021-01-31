/*
 *
 */
#define BER	1e-7   //1e-16 (0), 1e-6
#define absErrorBound	0.000001 //default 0.0001=2^{-12} (-13?), 0.000001=2^{-20}, 0.00001=2^{-16}, 0.001=2^{-10}, 0.01=2^{-7}
// #define absErrorBound_binary  20 //bitwise, SZ, equal to above
// #define relBoundRatio       0.01
// #define pw_relBoundRatio    0.01
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
*/
#define CT	10 
#define byte_or_bit         2 //1 byte, 2 bit
//#define data_num            8192 //pingpong
#define filename            "dataset/obs_info" //pingpong, k-means, "input", "testfloat_8_8_128", "testdouble_8_8_128", "testdouble_8_8_8_128", test, obs_info, num_plasma
#define suffix              ".txt" //k-means, ".txt"
#define output_suffix       "_output_" //k-means, "_output_", "_output_s_"
#define clusters            100 //k-means
//sz
#define bin_suffix          ".dat"
#define sz_suffix           ".sz"
#define zs_suffix           ".zs"
#define out_suffix           ".out"
#define sz_comp_cmd_prefix  "./sz -z -f -c sz.config -M ABS -A "
#define sz_comp_cmd_prefix_double  "./sz -z -d -c sz.config -M ABS -A "
#define sz_comp_cmd_suffix1 " -i "
#define sz_comp_cmd_suffix2 ".dat -1 "
#define sz_decomp_cmd_prefix  "./sz -x -f -s "
#define sz_decomp_cmd_prefix_double  "./sz -x -d -s "
#define sz_decomp_cmd_suffix ".dat.zs -1 "

void MPI_Bcast_bitwise_crc_hamming(double *buffer, int count, int root, int rank, int procs, float* compress_ratio, double* gosa, int* resend);
void MPI_Bcast_bitwise_mask_crc(double *buffer, int count, int root, int rank, int procs, float* compress_ratio, double* gosa, int* resend);
void MPI_Bcast_bitwise_crc(double *buffer, int count, int root, int rank, int procs, float* compress_ratio, double* gosa, int* resend);

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

float calCompressRatio_bitwise_float(float[], int);
float calCompressRatio_bitwise_double(double[], int);
float calCompressRatio_bitwise_double2(float[], int);

float calcCompressionRatio_sz_float(float[], int);
float calcCompressionRatio_nolossy_performance_float(float[], int);
float calcCompressionRatio_nolossy_area_float(float[], int);
float calcCompressionRatio_sz_double(double[], int);
float calcCompressionRatio_nolossy_performance_double(double[], int);
float calcCompressionRatio_nolossy_area_double(double[], int);

float calcCompressionRatio_himeno_ij_ik_jk(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
// MPI_Datatype myCompress_himeno(void*, int, int, int, int, int, int);
float calcCompressionRatio_himeno_sz(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
float calcCompressionRatio_himeno_nolossy_performance(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
float calcCompressionRatio_himeno_nolossy_area(float[MIMAX][MJMAX][MKMAX], int, int, int, int, int);
void getFloatBin(float, char[]);
void getDoubleBin(double,char[]);
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
void writetobinary_char(const char *, unsigned char*, int);
float* readfrombinary_float(const char*, int);
double* readfrombinary_double(const char*, int);
unsigned char* readfrombinary_char(const char *, int*);
float* readfrombinary_writetotxt_float(const char*, const char*, int);
double* readfrombinary_writetotxt_double(const char*, const char*, int);

void add_bit_to_bytes(unsigned char**, int*, int*, int);
void bit_set(unsigned char*, unsigned char, int);

int to_absErrorBound_binary(double absErrBound);

uint32_t do_crc32(unsigned char *data_bits, int bytes);
uint64_t get_random_int(uint64_t from, uint64_t to);

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
void bit_flip(unsigned char* bits, int bytes); //bit flip
int block_size(int data_bytes); //calculate block size (bytes)