# data-compression
This repo contains the work on MPI lossy floating-point data compression on the applications of Pingpong, Himeno, K-means, etc. 

## Usage: 
```shell
$ make
$ mpirun -np 2 ./pingpong
```

## Compile (one by one)
### Pingpong
```shell
$ mpicc pingpong.c dataCompression.c -o pingpong -lm -lz
```
### Himeno
```shell
$ mpicc himenoBMTxps.c dataCompression.c -o bmt -lm -lz
```
### K-means
```shell
$ mpicc k-means.c dataCompression.c -o k-means -lm -lz
```

## Execution 
### Pingpong
```shell
$ mpirun -np 2 ./pingpong
```
### Himeno
```shell
$ mpirun -np 2 ./bmt
```
### K-means
```shell
$ mpirun -np 2 ./k-means 2
```

## Source Files
### dataCompression.c
This file includes various floating-point data compression methods and other related functions for MPI communication.

### dataCompression.h
This file is the head file of dataCompression.c. The most important parameter to be used is CT.

| Parameter | Note | Value |
| --- | --- | --- |
| **CT** | Compression type |  0: no compress, <br> 1: byte-wise compression, <br> 2: no-lossy-performance compression, <br> 3: no-lossy-area compression, <br> 4: sz, <br> 5: bit-wise compression, <br> 6: bit-wise compression with no prediction, <br> 7 bitmask-based bit-wise compression, <br> 8 bit-wise compression with CRC |
| **absErrorBound** | Absolute error bound | Any |
| **clusters** | Number of clusters in K-means clustering | Integer (defaultly 100) |
| **filename** | file of raw data | dataset/ |
