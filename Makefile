CFLAGS = -Wall -O3
LDFLAGS = -lm -lfftw3
DEBUG =
D1 = 0
D2 = 0
SUB =

all: fftdiff interleave fftdiff_rect

fftdiff_rect: fftdiff.c
	$(CC) $(CFLAGS) $(LDFLAGS) fftdiff.c -o fftdiff_rect \
          -DPRESET -DWINDOW=win_rect -DTBUFSIZE=4096 -DPREBUFSIZE=1024 -DPOSTBUFSIZE=1024

clean:
	rm -f fftdiff interleave

test: fftdiff$(SUB) tdir/in1.raw tdir/in2.raw tdir/diff.raw
	$(DEBUG) ./fftdiff$(SUB) tdir/in1.raw $(D1) tdir/in2.raw $(D2) tdir/out.raw

tdir/diff.raw: tdir/in1.raw tdir/in2.raw
	./interleave tdir/in1.raw tdir/in2.raw tdir/diff.raw

%.raw: %.mp3
	mplayer -ao pcm -af volume=-3 -channels 1 -nowaveheader -aofile "$@" "$<"
