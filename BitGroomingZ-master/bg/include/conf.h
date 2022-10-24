/**
 *  @file conf.h
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief Header file for the conf.c.
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _Conf_H
#define _Conf_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdio.h>
#include "defines.h"

//conf.c
int BG_ReadConf(const char* bg_cfgFile);
int BG_LoadConf(const char* bg_cfgFile);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _Conf_H  ----- */

