/**
 *  @file   conf.c
 *  @author Sheng Di (sdi1@anl.gov or disheng222@gmail.com)
 *  @date   2015.
 *  @brief  Configuration loading functions for the BG library.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <math.h>
#include "string.h"
#include "iniparser.h"
#include "conf.h"
#include "bg.h"

/*-------------------------------------------------------------------------*/
/**
    @brief      It reads the configuration given in the configuration file.
    @return     integer         1 if successfull.

    This function reads the configuration given in the BG configuration
    file and sets other required parameters.

 **/
 
/*struct node_t *pool;
node *qqq;
node *qq;
int n_nodes = 0, qend;
unsigned long **code;
unsigned char *cout;
int n_inode;*/ 
 
/*-------------------------------------------------------------------------*/
/**
 * 
 * 
 * @return the status of loading conf. file: 1 (success) or 0 (error code);
 * */
int BG_ReadConf(const char* bg_cfgFile) {
    // Check access to BG configuration file and load dictionary
    //record the setting in confparams_cpr
    confparams_cpr = (bg_params*)malloc(sizeof(bg_params));    
    exe_params = (bg_exedata*)malloc(sizeof(bg_exedata));
    
    int x = 1;
    char sol_name[256];
    char *modeBuf;
    char *endianTypeString;
    dictionary *ini;
    char *par;

    char *y = (char*)&x;
	
    if(*y==1)
	sysEndianType = LITTLE_ENDIAN_SYSTEM;
    else //=0
	sysEndianType = BIG_ENDIAN_SYSTEM;
    
    if(bg_cfgFile == NULL)
    {
		dataEndianType = LITTLE_ENDIAN_DATA;
		confparams_cpr->sol_ID = BG;
				
		confparams_cpr->zlibMode = 1; //high speed mode
		
		return BG_SCES;
	}
    
    if (access(bg_cfgFile, F_OK) != 0)
    {
        printf("[BG] Configuration file NOT accessible.\n");
        return BG_NSCS;
    }
    
    ini = iniparser_load(bg_cfgFile);
    if (ini == NULL)
    {
        printf("[BG] Iniparser failed to parse the conf. file.\n");
        return BG_NSCS;
    }

	endianTypeString = iniparser_getstring(ini, "ENV:dataEndianType", "LITTLE_ENDIAN_DATA");
	if(strcmp(endianTypeString, "LITTLE_ENDIAN_DATA")==0)
		dataEndianType = LITTLE_ENDIAN_DATA;
	else if(strcmp(endianTypeString, "BIG_ENDIAN_DATA")==0)
		dataEndianType = BIG_ENDIAN_DATA;
	else
	{
		printf("Error: Wrong dataEndianType: please set it correctly in bg.config.\n");
		iniparser_freedict(ini);
		return BG_NSCS;
	}

	// Reading/setting detection parameters
	
	par = iniparser_getstring(ini, "ENV:sol_name", NULL);
	snprintf(sol_name, 256, "%s", par);
	
    if(strcmp(sol_name, "BG")==0)
		confparams_cpr->sol_ID = BG;
    else{
		printf("[BG] Error: wrong solution name (please check bg.config file), sol=%s\n", sol_name);
		iniparser_freedict(ini);
		return BG_NSCS;
    }
	
	modeBuf = iniparser_getstring(ini, "PARAMETER:zlibMode", NULL);
	if(modeBuf==NULL)
	{
		printf("[BG] Error: Null Zlib mode setting (please check bg.config file)\n");
		iniparser_freedict(ini);
		return BG_NSCS;					
	}		
	else if(strcmp(modeBuf, "Zlib_NO_COMPRESSION")==0)
		confparams_cpr->zlibMode = 0;
	else if(strcmp(modeBuf, "Zlib_BEST_SPEED")==0)
		confparams_cpr->zlibMode = 1;
	else if(strcmp(modeBuf, "Zlib_BEST_COMPRESSION")==0)
		confparams_cpr->zlibMode = 9;
	else if(strcmp(modeBuf, "Zlib_DEFAULT_COMPRESSION")==0)
		confparams_cpr->zlibMode = -1;
	else
	{
		printf("[BG] Error: Wrong zlib Mode /please check bg.config file)\n");
		return BG_NSCS;
	}

	modeBuf = iniparser_getstring(ini, "PARAMETER:bgMode", "BITGROOM");
	if(strcmp(modeBuf, "BITGROOM")==0)
		confparams_cpr->bgMode = BITGROOM;
	else if(strcmp(modeBuf, "BITSHAVE")==0)
		confparams_cpr->bgMode = BITSHAVE;
	else if(strcmp(modeBuf, "BITSET")==0)
		confparams_cpr->bgMode = BITSET;
	
	modeBuf = iniparser_getstring(ini, "PARAMETER:errorControlMode", BG_NSD);
	if(strcmp(modeBuf, "NSD")==0)
		confparams_cpr->errorControlMode = BG_NSD;
	else if(strcmp(modeBuf, "DSD")==0)
		confparams_cpr->errorControlMode = BG_DSD;
	
	confparams_cpr->NSD = (int)iniparser_getint(ini, "PARAMETER:NSD", 0);
	confparams_cpr->DSD = (int)iniparser_getint(ini, "PARAMETER:DSD", 0);

    iniparser_freedict(ini);
    return BG_SCES;
}

/*-------------------------------------------------------------------------*/
/**
    @brief      It reads and tests the configuration given.
    @return     integer         1 if successfull.

    This function reads the configuration file. Then test that the
    configuration parameters are correct (including directories).

 **/
/*-------------------------------------------------------------------------*/
int BG_LoadConf(const char* bg_cfgFile) {
    int res = BG_ReadConf(bg_cfgFile);
    if (res != BG_SCES)
    {
        printf("[BG] ERROR: Impossible to read configuration.\n");
        return BG_NSCS;
    }
    return BG_SCES;
}
