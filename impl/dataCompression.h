/*
 *
 */
#define absErrBound         0.0001
#define relBoundRatio       0.01
#define pw_relBoundRatio    0.01    
#define CT                  2 //compress type

float calcCompressionRatio_himeno_ij_ik_jk(float[MIMAX][MJMAX][MKMAX], int, int);
// MPI_Datatype myCompress_himeno(void*, int, int, int, int, int, int);

float calcCompressionRatio_himeno_nolossy_performance(float[MIMAX][MJMAX][MKMAX], int, int);
void getFloatBin(float, char[]);