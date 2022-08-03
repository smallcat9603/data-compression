#!/bin/sh

for datasize in 8192 16384 32768 65536 131072 262144
do
    for compression in bytewise bitwise bitmask bitnp bitop
    do
        echo "./mycompress_${compression}_double float_eq_${datasize}.txt" >> run.txt
        ./mycompress_${compression}_double float_eq_${datasize}.txt >> run.txt
    done 
done