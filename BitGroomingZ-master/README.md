BGZ: Bit Grooming Compressor (the Bit Grooming code extracted from Zender's nco, for the convenience of testing its compresibility)
=====
 (C) 2020 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
       See COPYRIGHT in top-level directory.

* Major Authors: Sheng Di 
* Supervisor: Franck Cappello 
* Other Contributors: 

## Installation

### Installation way 1:
* ./configure --prefix=[INSTALL_DIR]
* make
* make install

Then, you'll find all the executables in [INSTALL_DIR]/bin and .a and .so libraries in [INSTALL_DIR]/lib

## Testing Examples
--------------------------------------

Examples can be found in the [SZ_PACKAGE]/examples


## Compression
--------------
* ./testfloat_compress sz.config testfloat_8_8_128.dat 8192

## Decompression

* ./testfloat_decompress testdouble_8_8_128.dat.sz 8192

