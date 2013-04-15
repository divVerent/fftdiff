#!/bin/sh

a=$1
i=$2
b=$3

rawspec="-t raw -e signed -b 16 -c 1"

case "$a" in
	"--synth="*)
		synth=${a#*=}
		r=${synth%% *}
		s=${synth#* }
		sox $rawspec -r "$r" -n $rawspec synth.raw synth $s
		channels=synth.raw
		merge=
		;;
	*)
		c=`soxi -c "$a"`
		r=`soxi -r "$a"`
		k=0
		channels=
		while [ $k -lt $c ]; do
			k=$(($k+1))
			sox "$a" $rawspec "$k.raw" remix "$k"
			channels="$channels $k.raw"
		done
		if [ $c -gt 1 ]; then
			merge=--combine=merge
		fi
		;;
esac

geometry=
set --
for c in $channels; do
	if [ -z "$geometry" ]; then
		geometry=`./fftimgeq "$c" /dev/null /dev/null 2>&1 | tee /dev/stderr | awk '/^Total image size:/ { print $4 }'`
	fi
	convert "$i" -flip -transpose -geometry "$geometry"! -depth 8 -grayscale Rec601Luma GRAY:- | ./fftimgeq "$c" /dev/stdin "$c".out
	set -- "$@" $rawspec -r "$r" "$c".out
	out="$out $c.out"
done

sox $merge "$@" "$b"
