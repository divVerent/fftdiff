#!/bin/sh

a=$1
i=$2
b=$3

r=44100 # TODO detect this

sox "$a" -e signed -b 16 -c 1 1.raw
geometry=`./fftimgeq 1.raw /dev/null 2.raw 2>&1 | awk '/^Total image size:/ { print $4 }'`
convert "$i" -flip -transpose -geometry "$geometry"! -depth 8 GRAY:- |\
	./fftimgeq 1.raw /dev/stdin 2.raw
sox -t raw -e signed -b 16 -c 1 -r "$r" 2.raw "$b"
