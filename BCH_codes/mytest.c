/*
 * Generic binary BCH encoding/decoding library
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 as published by
 * the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 51
 * Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * This library is a fork of http://lxr.free-electrons.com/source/lib/bch.c
 *
 * Author: Alessandro Budroni
 *
 * Modifications: Creation of two functions to encode/decode 128-bit messages.
 *
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "bch_functions.h"

int main(){

    int len_msg = 16; //16; //bytes
    int correctable_errors = 4; //bytes
    int LOOP = 1000;
    
    time_t t;
    int i,j, count = LOOP, pEncodedLen, pDecodedLen, errors, ret, maxEncodedLen = len_msg + correctable_errors;
    unsigned char msg[len_msg], original_msg[len_msg];
    unsigned char pEncoded[maxEncodedLen], pDecoded[len_msg];
    
    srand((unsigned) time(&t));
    
    printf("Start testing.\n\n");
    
    printf("Test correcting 0-%d errors:...           ",maxEncodedLen-len_msg);
    for (j = 1; j<=LOOP; j++) {
    
        // Generate message
    
        for (i = 0; i< len_msg; i++) {
            msg[i] = rand() % 256;
        }
        
        memcpy(original_msg, msg, len_msg);
    
        // Generate the code
        
        if(len_msg == 16){
            GenerateBCH128( msg, len_msg, pEncoded, maxEncodedLen, &pEncodedLen);
        }
        else{
            GenerateBCH( msg, len_msg, pEncoded, maxEncodedLen, &pEncodedLen, correctable_errors);
        }     

        // printf("pEncoded: %d \n", pEncodedLen);
        
        // Generate errors
    
        errors = rand() % (maxEncodedLen-len_msg+1);
    
        for (i= 0; i<errors; i++){
            pEncoded[rand() % (pEncodedLen+1)] ^= ((unsigned char) 1 << ((rand() % 8)) & 0xFF);
        }
    
        // Decode and correct
    
        if(len_msg == 16){
            ValidateBCH128(pEncoded, pEncodedLen, pDecoded, maxEncodedLen, &pDecodedLen);
        }
        else{
            ValidateBCH(pEncoded, pEncodedLen, pDecoded, maxEncodedLen, &pDecodedLen, correctable_errors);
        }
    
        if(strncmp((char*)original_msg, (char*)pDecoded, len_msg) != 0){
            count--;
            printf("     - failed with %d errors test n. %d\n",errors,j-1);
        }
    }
    printf("\nfinished %d tests, %.1f%% passed!\n\n",j-1, 100.0*count/LOOP);
    
    count = LOOP;
    printf("Test detecting more than %d errors:...    \n", maxEncodedLen-len_msg);
    for (j = 1; j<=LOOP; j++) {
        
        // Generate message
        
        for (i = 0; i< len_msg; i++) {
            msg[i] = rand() % 256;
        }
        
        memcpy(original_msg, msg, len_msg);
        
        // Generate the code
        
        if(len_msg == 16){
            GenerateBCH128( msg, len_msg, pEncoded, maxEncodedLen, &pEncodedLen);
        }
        else{
            GenerateBCH( msg, len_msg, pEncoded, maxEncodedLen, &pEncodedLen, correctable_errors);
        }
        
        // Generate errors
        
        errors = rand() % 50;
        
        for (i= 0; i<errors; i++){
            pEncoded[rand() % (pEncodedLen+1)] ^= ((unsigned char) 1 << ((rand() % 8)) & 0xFF);
        }
        
        // Decode and correct
        
        if(len_msg == 16){
            ret = ValidateBCH128(pEncoded, pEncodedLen, pDecoded, maxEncodedLen, &pDecodedLen);
        }
        else{
            ValidateBCH(pEncoded, pEncodedLen, pDecoded, maxEncodedLen, &pDecodedLen, correctable_errors);
        } 
        
        if(strncmp((char*)original_msg, (char*)pDecoded, len_msg) != 0)
        {
            count--;
            printf("     - failed with %d errors, test n. %d\n",errors,j-1);
        }
    }
    printf("finished %d tests, %.1f%% passed!\n\n", j-1, 100.0*count/LOOP);
    printf("End testing.\n");
    
    return 0;
}



