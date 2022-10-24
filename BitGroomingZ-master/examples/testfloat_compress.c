/**
 *  @file test_compress.c
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief This is an example of using compression interface
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "bg.h"
#include "rw.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
	totalCost = 0;
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    int status = 0;
    size_t nbEle = 0;
    char oriFilePath[640], outputFilePath[640];
    char *cfgFile;
    
    if(argc < 3)
    {
		printf("Test case: testfloat_compress [config_file] [srcFilePath] [# elements]\n");
		printf("Example: testfloat_compress bg.config testfloat_8_8_128.dat 8192\n");
		exit(0);
    }
   
    cfgFile=argv[1];
    sprintf(oriFilePath, "%s", argv[2]);
    nbEle = atoi(argv[3]); //8
   
    printf("cfgFile=%s\n", cfgFile); 
    sprintf(outputFilePath, "%s.bg", oriFilePath);
   
	BG_LoadConf(cfgFile);
   
    float *data = readFloatData(oriFilePath, &nbEle, &status);
    if(status != BG_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
   
    size_t outSize = 1; 
    cost_start();
    unsigned char *bytes = BG_compress(BG_FLOAT, data, &outSize, nbEle);
    cost_end();
    printf("timecost=%f\n",totalCost); 
    writeByteData(bytes, outSize, outputFilePath, &status);
    if(status != BG_SCES)
    {
        printf("Error: data file %s cannot be written!\n", outputFilePath);
        exit(0);
    }

    printf("done\n");
    free(bytes); 
    free(data);
    
    return 0;
}
