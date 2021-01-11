#!/bin/sh

prog=k-means
procs=4

for CT in 8
do
    for AEB in 0.000001 0.00001 0.0001 0.001 0.01
    do
        for BER in 1e-16 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5
        do
            sh ./set-parameter.sh dataCompression.h $CT $AEB $BER
            mpicc $prog.c dataCompression.c -o $prog -lm -lz
            mpirun -np $procs ./$prog $procs
        done
    done 
done