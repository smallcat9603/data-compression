#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libdround.h>
#include <rw.h>
#include "libdround.h"

#include <time.h>

int main(int argc, char* argv[])
{
	char oriFilePath[640], cmpFile[640], outputFilePath[640];
 	int prec;

	double start_time_comp, end_time_comp, start_time_decomp, end_time_decomp;

  	if(argc < 3)
  	{
        printf("Test case: digitfloat_CompDecomp [data type] [precision] [srcFilePath]\n");
        printf("Example: digitfloat_CompDecomp -f 1 testdata/x86/testfloat_8_8_128.dat\n");
        exit(0);
  	}

  	char* type = argv[1];
  	prec = atoi(argv[2]);

  	sprintf(oriFilePath, "%s", argv[3]);
	sprintf(cmpFile, "%s.dg", oriFilePath);
	sprintf(outputFilePath, "%s.dg.out", oriFilePath);

	size_t nbEle = 0;
	int status = 0;
	float* data = readFloatData(oriFilePath, &nbEle, &status);
	//for(i=0;i<10;i++)
	//	printf("data[%d]=%f\n", i, data[i]);	
	//start compress: the entireu compression includes three steps: digit rounding + bit shuffle + zlib
	unsigned long outSize = 0;
	void* out = NULL;
	if(strcmp(type, "-f")==0)
	{
		start_time_comp = clock();
		unsigned char* compressBytes = dround_compress(DROUND_FLOAT, data, nbEle, prec, &outSize);
		end_time_comp = clock();
		writeByteData(compressBytes, outSize, cmpFile, &status);
		start_time_decomp = clock();
		out = dround_decompress(DROUND_FLOAT, compressBytes, nbEle, outSize);
		end_time_decomp = clock();
		writeFloatData_inBytes(out, nbEle, outputFilePath, &status);
	}
	else if(strcmp(type, "-d")==0)
	{
		start_time_comp = clock();
		unsigned char* compressBytes = dround_compress(DROUND_DOUBLE, data, nbEle, prec, &outSize);
		end_time_comp = clock();
		writeByteData(compressBytes, outSize, cmpFile, &status);
		start_time_decomp = clock();
		out = dround_decompress(DROUND_DOUBLE, compressBytes, nbEle, outSize);
		end_time_decomp = clock();
		writeDoubleData_inBytes(out, nbEle, outputFilePath, &status);
	}

	printf("Compression time: %f sec \n", (end_time_comp-start_time_comp)/CLOCKS_PER_SEC); 
	printf("Decompression time: %f sec \n", (end_time_decomp-start_time_decomp)/CLOCKS_PER_SEC); 
	
	//for(i=0;i<10;i++)
	//	printf("out[%d]=%f\n", i, out[i]);	
	// //free memory
	free(data);
	free(out);
}
