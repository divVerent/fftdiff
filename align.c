#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <err.h>

#include "types.h"

#define FALL 0.99

void getvol(char *inf, off_t ino, size_t len, double out[])
{
  FILE *f = fopen(inf, "rb");
  double cur = 0;
  int16_t icur;
  off_t i;

  if(!f)
    err(1, "fopen %s", inf);

  for(i = 0; i != ino; ++i)
  {
    fread(&icur, 1, sizeof(icur), f);
    cur = FALL * cur + (1 - FALL) * (double)icur * (double)icur;
  }
  for(i = 0; i != len; ++i)
  {
    fread(&icur, 1, sizeof(icur), f);
    cur = FALL * cur + (1 - FALL) * (double)icur * (double)icur;
    *out++ = cur;
  }

  fclose(f);
}

double squaresum(double a[], double b[], size_t n)
{
  double s = 0.0;
  while(n--)
    s += pow(*a++ - *b++, 2);
  return s;
}

void aligncheck(char *inf1, off_t ino1, char *inf2, off_t ino2, size_t len, size_t maxdelta)
{
  double *vol1 = malloc((len + maxdelta) * sizeof(double));
  double *vol2 = malloc((len + maxdelta) * sizeof(double));
  ssize_t delta, bestdelta = -2342;
  double bestscore = -1;

  if(!vol1 || !vol2)
    err(1, "malloc");

  getvol(inf1, ino1, len + maxdelta, vol1);
  getvol(inf2, ino2, len + maxdelta, vol2);

  for(delta = 0; delta <= maxdelta; ++delta)
  {
    double s = squaresum(vol1 + delta, vol2, len);
    if(bestscore < 0 || s < bestscore)
    {
      bestscore = s;
      bestdelta = delta;
    }
  }
  for(delta = 0; delta <= maxdelta; ++delta)
  {
    double s = squaresum(vol1, vol2 + delta, len);
    if(bestscore < 0 || s < bestscore)
    {
      bestscore = s;
      bestdelta = -delta;
    }
  }

  printf("Best delta: %ld (score: %f)\n", (long)bestdelta, bestscore);

  free(vol2);
  free(vol1);
}

int main(int argc, char **argv)
{
  if(argc != 7)
  {
    fprintf(stdout, "Usage: %s infile offset instrfile offset length maxdelta\n", *argv);
    fprintf(stdout, "  checks the alignment of the files\n");
    return 1;
  }
  aligncheck(argv[1], atol(argv[2]), argv[3], atol(argv[4]), atol(argv[5]), atol(argv[6]));
  return 0;
}
