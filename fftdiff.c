#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#include <err.h>

#include <fftw3.h>

#include "types.h"

#ifndef PRESET
#define TBUFSIZE 8192
#define PREBUFSIZE 3574
#define POSTBUFSIZE 3574
#define WINDOW win_blackman
#define CROSSFADE crf_triangle
#endif

#define MAXDEV 100

#define PBUFSIZE (PREBUFSIZE + POSTBUFSIZE)
#define BUFSIZE (TBUFSIZE - PBUFSIZE)
#define FTBUFSIZE (TBUFSIZE / 2 + 1) 

typedef struct
{
  FILE *data;
  wavedata buf[TBUFSIZE];
  int first;
} infile;

infile *infile_open(char *name, off_t pos)
{
  infile *a = malloc(sizeof(infile));
  if(!a) errx(1, "malloc: out of memory");
  if(!(a->data = fopen(name, "rb")))
    err(1, "open infile");
  fseek(a->data, pos * sizeof(wavedata), SEEK_SET);
  memset(a->buf, 0, sizeof(a->buf));
  a->first = 1;
  return a;
}

int infile_read(infile *f)
{
  size_t nread;

  if(f->first)
  {
    /* 1. read in a full buffer */
    nread = fread(f->buf, 1, sizeof(f->data), f->data);
    /* 2. never again */
    f->first = 0;
  }
  else
  {
    /* 1. shift the data to the left */
    memmove(f->buf, f->buf + BUFSIZE, PBUFSIZE * sizeof(f->buf[0]));
    /* 2. read in the missing data */
    nread = fread(f->buf + PBUFSIZE, 1, BUFSIZE * sizeof(f->buf[0]), f->data);
  }

  return nread;
}

windata windowing[TBUFSIZE];
windata crossfading[TBUFSIZE];
void win_hanning()
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
  {
    double a = 2*M_PI*i/(TBUFSIZE - 1);
    windowing[i] = (windata) ((MAX_WINDATA - 1.0) * (0.5 - 0.5*cos(a))) + 1;
  }
}
void win_blackman()
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
  {
    double a = 2*M_PI*i/(TBUFSIZE - 1);
    windowing[i] = (windata) ((MAX_WINDATA - 1.0) * (0.42 - 0.5*cos(a) + 0.08*cos(2*a))) + 1;
  }
}
void win_rect()
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
    windowing[i] = MAX_WINDATA;
}
void win_check()
{
  size_t i;
  windata min = MAX_WINDATA;
  for(i = 0; i != TBUFSIZE; ++i)
  {
    if(i >= PREBUFSIZE && i < PREBUFSIZE + BUFSIZE)
      if(min > windowing[i])
        min = windowing[i];
  }
  printf("windowing quality: %u\n", (unsigned)min);
}
void crf_none()
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
  {
    crossfading[i] =
      (i < PREBUFSIZE || i >= PREBUFSIZE + BUFSIZE) ? 0 : 65535;
  }
}
void crf_triangle()
{
  size_t i;
  const size_t left = PREBUFSIZE - BUFSIZE/2;
  const size_t middle = PREBUFSIZE + BUFSIZE/2;
  const size_t right = PREBUFSIZE + BUFSIZE + BUFSIZE/2;
  for(i = 0; i != TBUFSIZE; ++i)
  {
    crossfading[i] =
      (i <= left || i >= right) ? 0 :
      (i <= middle) ? (i  - left) * (winwindata)MAX_WINDATA / (middle  - left) :
      (i >  middle) ? (right - i) * (winwindata)MAX_WINDATA / (right - middle) :
      -2342;
  }
}
void crf_fix()
{
  size_t i, r;
  for(r = 0; r != BUFSIZE; ++r)
  {
    winwindata s = 0;
    windata S = MAX_WINDATA;
    for(i = r; i < TBUFSIZE; i += BUFSIZE)
    {
      s += crossfading[i];
    }
    if(s == 0)
      errx(1, "Your crossfading window is broken, fix it! at winpos %u", (unsigned)r);
    if(s > S + MAXDEV || s + MAXDEV < S)
      warnx("Your window seems to be broken, keep out of reach of children! %u != %u at %u", (unsigned)s, (unsigned)S, (unsigned)r);
    for(i = r; i < TBUFSIZE && s; i += BUFSIZE)
    {
      wavwindata sold = s;
      s -= crossfading[i];
      crossfading[i] = (crossfading[i] * (winwindata)S) / sold;
      S -= crossfading[i];
    }
    if(s != 0 || S != 0)
      errx(1, "I am broken, please fix me. Now. %u vs. %u at winpos %d", (unsigned)s, (unsigned)S, r);
  }
}


void buf_preprocess(wavedata buf[], wavedata obuf[])
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
    obuf[i] = (((wavwindata) buf[i]) * windowing[i]) / MAX_WINDATA;
}
void buf_postprocess(wavedata buf[])
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
    buf[i] = (((wavwindata) buf[i]) * MAX_WINDATA) / windowing[i];
}
void buf_crossfade(wavedata inbuf[], wavedata outbuf[])
{
  size_t i;
  /* 1. Shift everything by BUFSIZE to the left. */
  memmove(outbuf, outbuf + BUFSIZE, PBUFSIZE * sizeof(*outbuf));
  memset(outbuf + PBUFSIZE, 0, BUFSIZE * sizeof(*outbuf));
  /* 2. write the data to the new buffer */
  for(i = 0; i != TBUFSIZE; ++i)
    outbuf[i] += inbuf[i] * (wavwindata)crossfading[i] / MAX_WINDATA;
}

typedef struct
{
  double *buf1, *buf2, *obuf;
  complex *fbuf1, *fbuf2, *fobuf;
  fftw_plan plan1, plan2, oplan;
} fftdata;

void fft_init(fftdata *d)
{
  d->buf1 = fftw_malloc(TBUFSIZE * sizeof(*d->buf1));
  d->buf2 = fftw_malloc(TBUFSIZE * sizeof(*d->buf2));
  d->obuf = fftw_malloc(TBUFSIZE * sizeof(*d->obuf));
  d->fbuf1 = fftw_malloc(FTBUFSIZE * sizeof(*d->fbuf1));
  d->fbuf2 = fftw_malloc(FTBUFSIZE * sizeof(*d->fbuf2));
  d->fobuf = fftw_malloc(FTBUFSIZE * sizeof(*d->fobuf));
  fputs("Preparing FFTs.", stderr);
  d->plan1 = fftw_plan_dft_r2c_1d(TBUFSIZE, d->buf1, d->fbuf1, 0);
  fputc('.', stderr);
  d->plan2 = fftw_plan_dft_r2c_1d(TBUFSIZE, d->buf2, d->fbuf2, 0);
  fputc('.', stderr);
  d->oplan = fftw_plan_dft_c2r_1d(TBUFSIZE, d->fobuf, d->obuf, 0);
  fputc('\n', stderr);
}

void fft_free(fftdata *d)
{
  fftw_destroy_plan(d->oplan);
  fftw_destroy_plan(d->plan2);
  fftw_destroy_plan(d->plan1);
  fftw_free(d->fobuf);
  fftw_free(d->fbuf2);
  fftw_free(d->fbuf1);
  fftw_free(d->obuf);
  fftw_free(d->buf2);
  fftw_free(d->buf1);
}

void d2w(double in[], wavedata out[])
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
    out[i] = in[i] / TBUFSIZE;
}

void w2d(wavedata in[], double out[])
{
  size_t i;
  for(i = 0; i != TBUFSIZE; ++i)
    out[i] = in[i];
}

void buf_convolve(fftdata *f, wavedata buf1[], wavedata buf2[], wavedata obuf[])
{
  size_t i;
  complex *b1 = f->fbuf1, *b2 = f->fbuf2, *ob = f->fobuf;
  w2d(buf1, f->buf1); fftw_execute(f->plan1);
  w2d(buf2, f->buf2); fftw_execute(f->plan2);
  for(i = 0; i != 20; ++i)
    ob[i] = b1[i];
  for(i = 0; i != FTBUFSIZE; ++i)
  {
    double a1 = cabs(b1[i]);
    double a = a1 - cabs(b2[i]) * 1.1;
    if(a <= 0)
      ob[i] = 0;
    else
      ob[i] = b1[i] * a / a1;
  }
  fftw_execute(f->oplan); d2w(f->obuf, obuf);
}
  
void fftdiff(char *in1, off_t pos1, char *in2, off_t pos2, char *out)
{
  fftdata f;
  wavedata buf[TBUFSIZE];
  wavedata ib1[TBUFSIZE];
  wavedata ib2[TBUFSIZE];
  wavedata obuf[TBUFSIZE] = {0};
  infile *if1 = infile_open(in1, pos1);
  infile *if2 = infile_open(in2, pos2);
  FILE *of = fopen(out, "wb");
  if(!of) err(1, "open output file");

  fft_init(&f);
  for(;;)
  {
    if(!infile_read(if1))
      break;
    if(!infile_read(if2))
      break;
    buf_preprocess(if1->buf, ib1);
    buf_preprocess(if2->buf, ib2);
    buf_convolve(&f, ib1, ib2, buf);
    buf_postprocess(buf);
    buf_crossfade(buf, obuf);
    if(!fwrite(obuf, 1, BUFSIZE * sizeof(obuf[0]), of))
      err(1, "writing");
    putc('.', stderr);
  }
  fft_free(&f);
  
  putc('\n', stderr);
}

int main(int argc, char **argv)
{
  if(argc != 6)
  {
    fprintf(stdout, "Usage: %s infile offset instrfile offset outfile\n", *argv);
    fprintf(stdout, "  will remove the instruments from instrfile from infile and write the result to\n");
    fprintf(stdout, "  outfile. File format: RAW, 16bit mono in host byte order.\n");
    return 1;
  }
  WINDOW();
  CROSSFADE();
  crf_fix();
  win_check();
  fftdiff(argv[1], atol(argv[2]), argv[3], atol(argv[4]), argv[5]);
  return 0;
}
