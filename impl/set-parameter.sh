#!/bin/sh
fn=$1
ct=$2
aeb=$3
ber=$4

sed -i "s/#define CT[[:space:]][[:space:]]*[^ \t][^ \t]*/#define CT\t$ct/" $fn
sed -i "s/#define absErrorBound[[:space:]][[:space:]]*[^ \t][^ \t]*/#define absErrorBound\t$aeb/" $fn
sed -i "s/#define BER[[:space:]][[:space:]]*[^ \t][^ \t]*/#define BER\t$ber/" $fn