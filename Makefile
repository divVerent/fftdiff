CFLAGS = -Wall -O3
LDFLAGS = -lm -lfftw3
DEBUG =

all: fftdiff interleave fftdiff_rect

fftdiff_rect: fftdiff.c
	$(CC) $(CFLAGS) $(LDFLAGS) fftdiff.c -o fftdiff_rect \
          -DPRESET -DWINDOW=win_rect -DTBUFSIZE=4096 -DPREBUFSIZE=1024 -DPOSTBUFSIZE=1024

clean:
	rm -f fftdiff interleave

test: fftdiff
	$(DEBUG) ./fftdiff test/in1.raw 0 test/in2.raw 0 test/out.raw
