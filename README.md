# data-compression
This repo contains the work on MPI lossy floating-point data compression on the applications of Pingpong, Himeno, K-means, FFT, MM, LU, etc. 

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

### MM
```shell
$ mpicc mm.c dataCompression.c -o mm -lm -lz
```

### LU
```shell
$ mpicc lu.c dataCompression.c -o lu -lm -lz
```

## Execution 
### Pingpong
```shell
$ mpirun -np 2 ./pingpong [CT]
```
### Himeno
```shell
$ mpirun -np 2 ./bmt [CT]
```
### K-means
```shell
$ mpirun -np 2 ./k-means [CT]
```

### MM
```shell
$ mpirun -np 4 ./mm testdata/mat_1000_1000_a.txt testdata/mat_1000_1000_b.txt [CT]
```

### LU
```shell
$ mpirun -np 4 ./lu 16 [CT]
```

## Source Files
### dataCompression.c
This file includes various floating-point data compression methods and other related functions for MPI communication.

### dataCompression.h
This file is the head file of dataCompression.c.

### set-parameter.sh
This file is used to set the parameters of **CT**, **absErrorBound** and **BER** in dataCompression.h.
#### Usage
```shell
$ sh ./set-parameter.sh dataCompression.h 8 0.01 1e-6 // CT = 8, absErrorBound = 0.01, BER = 1e-6
```

| Parameter | Note | Value |
| --- | --- | --- |
| **CT** | Compression type |  0: no compress, <br> 1: byte-wise compression, <br> 2: no-lossy-performance compression, <br> 3: no-lossy-area compression, <br> 4: sz, <br> 5: bit-wise compression, <br> 6: bit-wise compression with no prediction, <br> 7 bitmask-based bit-wise compression, <br> 8 bit-wise compression with CRC, <br> 9 bitmask-based bit-wise compression with CRC, <br> 10 bit-wise compression with CRC Hamming, <br> 11 bit-wise compression with only prediction |
| **absErrorBound** | Absolute error bound | Any |
| **clusters** | Number of clusters in K-means clustering | Integer (defaultly 100) |
| **filename** | file of raw data | dataset/ |
