# data-compression
This repo contains the work on MPI lossy floating-point data compression on the applications of Pingpong, Himeno, K-means, etc. 

## Compile
* make (by Makefile) or compile one by one
* (Pingpong) mpicc pingpong.c dataCompression.c -o pingpong -lm
* (Himeno) mpicc himenoBMTxps.c dataCompression.c -o bmt -lm
* (K-means) mpicc k-means.c dataCompression.c -o k-means -lm

## Execution 
* (Pingpong) mpirun -np 2 ./pingpong
* (Himeno) mpirun -np 2 ./bmt
* (K-means) mpirun -np 2 ./k-means 2

