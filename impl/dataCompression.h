/*
 *
 */
#define absErrBound         0.0001 //2^{-12}
#define absErrBound_binary  12 //equal to above
#define relBoundRatio       0.01
#define pw_relBoundRatio    0.01    
#define CT                  4 //compress type
#define byte_or_bit         1 //1 byte, 2 bit

float calcCompressionRatio_himeno_ij_ik_jk(float[MIMAX][MJMAX][MKMAX], int, int);
// MPI_Datatype myCompress_himeno(void*, int, int, int, int, int, int);
float calcCompressionRatio_himeno_sz(float[MIMAX][MJMAX][MKMAX], int, int);
float calcCompressionRatio_himeno_nolossy_performance(float[MIMAX][MJMAX][MKMAX], int, int);
float calcCompressionRatio_himeno_nolossy_area(float[MIMAX][MJMAX][MKMAX], int, int);
void getFloatBin(float, char[]);
float* readFileFloat(char[]);