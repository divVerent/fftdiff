#!/bin/sh

a=$1
b=$2

set -ex

ofs=50000
len=32768
maxdelta=16384

sox "$a" -e signed -b 16 -c 1 1.raw
sox "$b" -e signed -b 16 -c 1 2.raw
r=`./align "$a" "$ofs" "$b" "$ofs" "$len" "$maxdelta" | tee /dev/stderr`
delta=`echo "$r" | cut -d ' ' -f 3`
case "$delta" in
	-*)
		./fftdiff 1.raw "0" 2.raw "${delta#-}" 3.raw
		;;
	*)
		./fftdiff 1.raw "$delta" 2.raw "0" 3.raw
		;;
esac
