CFLAGS = -Wall -ggdb3
LDFLAGS = -lm -lfftw3
DEBUG =

all: fftdiff interleave

test: fftdiff
	$(DEBUG) ./fftdiff
